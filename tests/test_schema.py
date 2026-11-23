import pyarrow as pa
import yaml
from pathlib import Path

from hk_sunu_sim.schema import TABLES

ROOT = Path(__file__).resolve().parents[1]


def test_schemas_have_unique_columns_and_batch():
    for name, schema in TABLES.items():
        names = schema.names
        assert len(names) == len(set(names)), name
        assert names[0] == "batch", name


def test_empty_tables_build():
    for schema in TABLES.values():
        pa.Table.from_pylist([], schema=schema)


def test_example_configs_load():
    for path in (ROOT / "config").glob("*.yaml"):
        cfg = yaml.safe_load(path.read_text())
        assert "production" in cfg and "seed" in cfg
