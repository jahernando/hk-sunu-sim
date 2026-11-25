#!/bin/bash -l
#SBATCH --job-name=genie-gevgen
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=%x-%A_%a.out

# GENIE event generation on water with the community container
# (docker://nuisancemc/genie:3.04.00, GENIE 3.04.00 + ROOT 6.30.04 + Pythia6).
# Each array task writes <out>/<tag>_<task>.ghep.root and .gst.root.
#
#   sbatch --array=0-9 slurm/genie_gevgen.sbatch
#   NU=-12 N=5000 TUNE=AR23_20i_00_000 sbatch --array=0-3 slurm/genie_gevgen.sbatch
#
# Environment (defaults):
#   PRODUCTION  atm_sk4_v0           output subdir under $LUSTRE/hk-sunu-sim
#   NU          14                   neutrino pdg
#   N           10000                events per task
#   TUNE        G18_10a_02_11b       one of the splines shipped in the image
#   EMIN,EMAX   0.05,10              GeV
#   FLUX        1/x                  gevgen -f: function of E [GeV] or file.root,hist
#   SEED_BASE   1                    seed = SEED_BASE*1000000 + task index
#   SIF         $STORE/software/hk/containers/genie_3.04.00.sif
#
# Water by mass: 16O 0.888, 1H 0.112. Splines in the image are the "min" set
# (nu_e, nu_mu and antineutrinos on 1H, 2H, 12C, 16O, 40Ar, 56Fe; 10 MeV-1 TeV);
# nu_tau needs its own splines (gmkspl).

module load cesga/2025 apptainer/1.4.2

PRODUCTION=${PRODUCTION:-atm_sk4_v0}
NU=${NU:-14}
N=${N:-10000}
TUNE=${TUNE:-G18_10a_02_11b}
EMIN=${EMIN:-0.05}
EMAX=${EMAX:-10}
FLUX=${FLUX:-1/x}
SEED_BASE=${SEED_BASE:-1}
SIF=${SIF:-$STORE/software/hk/containers/genie_3.04.00.sif}

TASK=${SLURM_ARRAY_TASK_ID:-0}
SEED=$((SEED_BASE * 1000000 + TASK))
TAG=genie_${TUNE}_nu${NU}
OUT=$LUSTRE/hk-sunu-sim/$PRODUCTION/genie
mkdir -p $OUT
WORK=${TMPDIR:-/tmp}/genie_$TASK          # gevgen writes input-flux.root in cwd: one dir per task
mkdir -p $WORK && cd $WORK

echo "host $(hostname)  task $TASK  seed $SEED  nu $NU  N $N  tune $TUNE  E $EMIN-$EMAX  flux $FLUX"

apptainer exec -B $LUSTRE,$WORK $SIF bash -c "
  XS=/opt/genie_xsec/3_04_00/$TUNE/gxspl-min.xml.gz
  gevgen -n $N -p $NU -t '1000080160[0.888],1000010010[0.112]' -e $EMIN,$EMAX -f '$FLUX' \
         --tune $TUNE --cross-sections \$XS --seed $SEED -r $TASK \
         -o ${TAG}_${TASK}.ghep.root > ${TAG}_${TASK}.gevgen.log 2>&1 || exit 1
  gntpc -i ${TAG}_${TASK}.ghep.root -f gst --tune $TUNE \
        -o ${TAG}_${TASK}.gst.root > ${TAG}_${TASK}.gntpc.log 2>&1 || exit 2
"
RC=$?
mv ${TAG}_${TASK}.* $OUT/
[ $RC -eq 0 ] || { echo "failed ($RC)"; tail -20 $OUT/${TAG}_${TASK}.gevgen.log; exit 1; }
ls -lh $OUT/${TAG}_${TASK}.*
