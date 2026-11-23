// wcprop: propagate primary particles in an infinite volume of water.
//
// Instead of tracking optical photons, the mean Cherenkov photon yield of every
// charged step is computed analytically (constant refractive index, 300-600 nm)
// and accumulated in 10 ns time bins. This "light" table is the input of the
// parametrised detector response. Neutron captures and the fate of stopping
// muons and pions are recorded as truth tables.
//
// Usage:
//   wcprop <primaries.csv> <output_prefix> [physics_list] [seed] [tmax_ns]
//
// primaries.csv: event_id,pdg,px,py,pz   (MeV/c; vertex at origin, t = 0)
// Outputs (CSV; cm, ns, MeV):
//   <prefix>_light.csv  event_id,tbin,n,sx,sy,sz,sux,suy,suz,st
//                       sums weighted by photon number n (x = step midpoint,
//                       u = direction at step start); mean = s*/n
//   <prefix>_ncap.csv   event_id,t,x,y,z,target_z,target_a,e_gamma
//   <prefix>_ends.csv   event_id,track_id,parent_id,pdg,process,fate,t,x,y,z,
//                       e_charged,n_neutron,n_gamma,e_gamma,n_nu,target_z
//                       one row per stopped/killed mu+-, pi+-; t is the time of
//                       the end products (decay or capture), fate is
//                       decay | capture | absorption | inelastic | other.
//                       n_gamma/e_gamma include muonic X-rays (~0.1-0.2 MeV)

#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserStackingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4PhysListFactory.hh"
#include "G4ThermalNeutrons.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cmath>

namespace {

struct Particle { int pdg; double px, py, pz; };
std::vector<long> gEventIds;
std::vector<std::vector<Particle>> gEvents;

std::ofstream gLight, gNcap, gEnds;
double gTmax = 2.0e6 * ns;
bool gThermal = false;

// Cherenkov constants: dN/dx = 2 pi alpha / (hc) * dE * sin^2(theta)
//   2 pi alpha / (hc) = 369.81 / (eV cm); dE = 1239.84 (1/300 - 1/600) eV
const double kRefIndex = 1.333;
const double kPhotonsPerCm = 369.81 * 1239.84 * (1.0 / 300.0 - 1.0 / 600.0);
const double kTimeBin = 10.0 * ns;

struct LightBin { double n = 0, sx = 0, sy = 0, sz = 0, sux = 0, suy = 0, suz = 0, st = 0; };
std::map<long, LightBin> gBins;

// End records are buffered during the event: the time of the end products is
// only final when they start being tracked (hadronic at-rest processes shift it
// after the step), so it is taken from the tracking action.
struct EndRecord { int track_id, parent_id, pdg, nn, ng, nnu, target_z; std::string process, fate;
                   double x, y, z, ech, eg; };
std::vector<EndRecord> gEndRecords;
std::map<int, double> gChildTime;            // ending track id -> earliest start time of its end products
std::map<const G4Track*, int> gEndChildren;  // end products (not delta rays) -> ending track id
long gCurrentEvent = -1;

void ReadPrimaries(const std::string& path) {
  std::ifstream in(path);
  if (!in) G4Exception("ReadPrimaries", "wcprop001", FatalException, ("cannot open " + path).c_str());
  std::string line;
  std::getline(in, line);  // header
  std::map<long, size_t> index;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string tok;
    std::vector<std::string> f;
    while (std::getline(ss, tok, ',')) f.push_back(tok);
    if (f.size() < 5) continue;
    long id = std::stol(f[0]);
    if (!index.count(id)) {
      index[id] = gEvents.size();
      gEventIds.push_back(id);
      gEvents.emplace_back();
    }
    gEvents[index[id]].push_back({std::stoi(f[1]), std::stod(f[2]), std::stod(f[3]), std::stod(f[4])});
  }
}

class Detector : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    G4Material* water;
    if (gThermal) {
      // Thermal scattering S(alpha, beta) for H bound in water needs this element name
      auto* hw = new G4Element("TS_H_of_Water", "H_WATER", 1.0, 1.0079 * g / mole);
      auto* o = G4NistManager::Instance()->FindOrBuildElement("O");
      water = new G4Material("Water_TS", 1.0 * g / cm3, 2, kStateLiquid, 293.15 * kelvin);
      water->AddElement(hw, 2);
      water->AddElement(o, 1);
    } else {
      water = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    }
    const double half = 100.0 * m;  // effectively infinite for MeV-GeV events
    auto* box = new G4Box("World", half, half, half);
    auto* lv = new G4LogicalVolume(box, water, "World");
    return new G4PVPlacement(nullptr, {}, lv, "World", nullptr, false, 0);
  }
};

class Primaries : public G4VUserPrimaryGeneratorAction {
 public:
  void GeneratePrimaries(G4Event* ev) override {
    size_t i = ev->GetEventID();
    gCurrentEvent = gEventIds[i];
    auto* vtx = new G4PrimaryVertex(G4ThreeVector(), 0.0);
    auto* table = G4ParticleTable::GetParticleTable();
    for (const auto& p : gEvents[i]) {
      G4ParticleDefinition* def = nullptr;
      if (p.pdg > 1000000000) {
        int z = (p.pdg / 10000) % 1000, a = (p.pdg / 10) % 1000;
        def = G4IonTable::GetIonTable()->GetIon(z, a, 0.0);
      } else {
        def = table->FindParticle(p.pdg);
      }
      if (!def) {
        G4cerr << "wcprop: skipping unknown pdg " << p.pdg << " in event " << gCurrentEvent << G4endl;
        continue;
      }
      auto* pp = new G4PrimaryParticle(def);
      pp->SetMomentum(p.px * MeV, p.py * MeV, p.pz * MeV);
      vtx->SetPrimary(pp);
    }
    ev->AddPrimaryVertex(vtx);
  }
};

class Stacking : public G4UserStackingAction {
 public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* t) override {
    return t->GetGlobalTime() > gTmax ? fKill : fUrgent;
  }
};

class Stepping : public G4UserSteppingAction {
 public:
  void UserSteppingAction(const G4Step* step) override {
    const G4Track* track = step->GetTrack();
    const auto* def = track->GetDefinition();
    const auto* pre = step->GetPreStepPoint();
    const auto* post = step->GetPostStepPoint();

    // Cherenkov yield of charged steps
    double q = def->GetPDGCharge() / eplus;
    if (q != 0.0 && step->GetStepLength() > 0) {
      double beta = 0.5 * (pre->GetBeta() + post->GetBeta());
      double s2 = 1.0 - 1.0 / (beta * beta * kRefIndex * kRefIndex);
      if (s2 > 0) {
        double n = kPhotonsPerCm * q * q * s2 * step->GetStepLength() / cm;
        double t = pre->GetGlobalTime();
        auto& b = gBins[static_cast<long>(t / kTimeBin)];
        G4ThreeVector x = 0.5 * (pre->GetPosition() + post->GetPosition()) / cm;
        G4ThreeVector u = pre->GetMomentumDirection();
        b.n += n;
        b.sx += n * x.x(); b.sy += n * x.y(); b.sz += n * x.z();
        b.sux += n * u.x(); b.suy += n * u.y(); b.suz += n * u.z();
        b.st += n * t / ns;
      }
    }

    const G4VProcess* proc = post->GetProcessDefinedStep();
    std::string pname = proc ? proc->GetProcessName() : "none";
    const auto* secs = step->GetSecondaryInCurrentStep();

    // Neutron captures
    if (def->GetPDGEncoding() == 2112 && pname == "nCapture") {
      double eg = 0; int tz = 0, ta = 0;
      for (const auto* s : *secs) {
        if (s->GetDefinition()->GetPDGEncoding() == 22) eg += s->GetKineticEnergy();
        else if (s->GetDefinition()->GetAtomicNumber() > 0) {
          tz = s->GetDefinition()->GetAtomicNumber();
          ta = s->GetDefinition()->GetAtomicMass() - 1;  // target = residual - captured n
        }
      }
      G4ThreeVector x = post->GetPosition() / cm;
      gNcap << gCurrentEvent << ',' << post->GetGlobalTime() / ns << ',' << x.x() << ',' << x.y() << ','
            << x.z() << ',' << tz << ',' << ta << ',' << eg / MeV << '\n';
    }

    // Fate of muons and charged pions
    int apdg = std::abs(def->GetPDGEncoding());
    if ((apdg == 13 || apdg == 211) && track->GetTrackStatus() == fStopAndKill) {
      double ech = 0, eg = 0; int nn = 0, ng = 0, nnu = 0;
      double tsec = post->GetGlobalTime();
      bool first = true;
      for (const auto* s : *secs) {
        const auto* sd = s->GetDefinition();
        int sp = sd->GetPDGEncoding();
        if (first || s->GetGlobalTime() < tsec) { tsec = s->GetGlobalTime(); first = false; }
        int asp = std::abs(sp);
        if (asp == 12 || asp == 14 || asp == 16) ++nnu;
        else if (sp == 2112) ++nn;
        else if (sp == 22) { ++ng; eg += s->GetKineticEnergy(); }
        else if (sd->GetPDGCharge() != 0 && sd->GetAtomicNumber() == 0) ech += s->GetKineticEnergy();
      }
      // The mu- at-rest process covers both decay in orbit and nuclear capture:
      // tell them apart by the neutrinos (2 in a decay, 1 in a capture).
      std::string fate = "other";
      if (pname.find("Inelastic") != std::string::npos) fate = "inelastic";
      else if (apdg == 13) fate = nnu >= 2 ? "decay" : (nnu == 1 ? "capture" : "other");
      else if (apdg == 211) fate = nnu >= 1 ? "decay" : "absorption";
      // Nucleus that took the particle (at-rest hadronic processes): needed to
      // correct the mu- atom assignment in water (see README of wcprop)
      int tz = 0;
      if (auto* hp = dynamic_cast<const G4HadronicProcess*>(proc))
        if (const G4Nucleus* nuc = hp->GetTargetNucleus()) tz = nuc->GetZ_asInt();
      G4ThreeVector x = post->GetPosition() / cm;
      gEndRecords.push_back({track->GetTrackID(), track->GetParentID(), def->GetPDGEncoding(), nn, ng, nnu, tz,
                             pname, fate, x.x(), x.y(), x.z(), ech / MeV, eg / MeV});
      gChildTime[track->GetTrackID()] = tsec;  // fallback, replaced by the tracking action
      // Timing: muonic X-rays and Auger electrons come out at stop time, before
      // the decay or capture. Use the neutrinos when there are any (decays,
      // mu- capture); otherwise (pi- absorption) all end products.
      for (const auto* sec : *secs) {
        int a = std::abs(sec->GetDefinition()->GetPDGEncoding());
        bool nu = (a == 12 || a == 14 || a == 16);
        if (nnu == 0 || nu) gEndChildren[sec] = track->GetTrackID();
      }
    }
  }
};

class Tracking : public G4UserTrackingAction {
 public:
  void PreUserTrackingAction(const G4Track* t) override {
    auto c = gEndChildren.find(t);
    if (c == gEndChildren.end()) return;
    int id = c->second;
    gEndChildren.erase(c);
    // first end product tracked replaces the fallback; keep the earliest afterwards
    if (!gTimeSet[id] || t->GetGlobalTime() < gChildTime[id]) {
      gChildTime[id] = t->GetGlobalTime();
      gTimeSet[id] = true;
    }
  }
  static std::map<int, bool> gTimeSet;
};
std::map<int, bool> Tracking::gTimeSet;

class Events : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event*) override {
    gBins.clear(); gEndRecords.clear(); gChildTime.clear(); gEndChildren.clear(); Tracking::gTimeSet.clear();
  }
  void EndOfEventAction(const G4Event*) override {
    for (const auto& r : gEndRecords)
      gEnds << gCurrentEvent << ',' << r.track_id << ',' << r.parent_id << ',' << r.pdg << ',' << r.process
            << ',' << r.fate << ',' << gChildTime[r.track_id] / ns << ',' << r.x << ',' << r.y << ',' << r.z
            << ',' << r.ech << ',' << r.nn << ',' << r.ng << ',' << r.eg << ',' << r.nnu << ',' << r.target_z << '\n';
    for (const auto& [k, b] : gBins)
      gLight << gCurrentEvent << ',' << k << ',' << b.n << ',' << b.sx << ',' << b.sy << ',' << b.sz << ','
             << b.sux << ',' << b.suy << ',' << b.suz << ',' << b.st << '\n';
  }
};

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    G4cerr << "usage: wcprop <primaries.csv> <output_prefix> [physics_list] [seed] [tmax_ns]" << G4endl;
    return 1;
  }
  std::string input = argv[1], prefix = argv[2];
  std::string plname = argc > 3 ? argv[3] : "QGSP_BIC_HP";
  long seed = argc > 4 ? std::stol(argv[4]) : 1;
  if (argc > 5) gTmax = std::stod(argv[5]) * ns;
  gThermal = true;  // always use S(alpha,beta) for H in water

  G4Random::setTheSeed(seed);
  ReadPrimaries(input);

  gLight.open(prefix + "_light.csv");
  gNcap.open(prefix + "_ncap.csv");
  gEnds.open(prefix + "_ends.csv");
  gLight << "event_id,tbin,n,sx,sy,sz,sux,suy,suz,st\n";
  gNcap << "event_id,t,x,y,z,target_z,target_a,e_gamma\n";
  gEnds << "event_id,track_id,parent_id,pdg,process,fate,t,x,y,z,e_charged,n_neutron,n_gamma,e_gamma,n_nu,target_z\n";

  auto* rm = new G4RunManager;
  rm->SetUserInitialization(new Detector);
  G4PhysListFactory factory;
  factory.SetVerbose(0);
  auto* pl = factory.GetReferencePhysList(plname);
  pl->RegisterPhysics(new G4ThermalNeutrons());  // radioactive decay is already in the _HP lists
  pl->SetVerboseLevel(0);
  rm->SetUserInitialization(pl);
  rm->SetUserAction(new Primaries);
  rm->SetUserAction(new Stacking);
  rm->SetUserAction(new Stepping);
  rm->SetUserAction(new Tracking);
  rm->SetUserAction(new Events);
  rm->Initialize();
  rm->BeamOn(static_cast<int>(gEvents.size()));

  gLight.close(); gNcap.close(); gEnds.close();
  delete rm;
  return 0;
}
