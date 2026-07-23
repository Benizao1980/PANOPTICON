from pathlib import Path
import importlib.util
import sys
import pandas as pd


def load_module():
    path = Path(__file__).parents[1] / "scripts" / "pangenome_io.py"
    spec = importlib.util.spec_from_file_location("pangenome_io", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_rtab_with_empty_csv(tmp_path):
    module = load_module()
    pd.DataFrame(
        {
            "Gene": ["famA", "famB", "famC"],
            "g1": [1, 1, 0],
            "g2": [1, 0, 1],
            "g3": [1, 0, 0],
        }
    ).to_csv(tmp_path / "gene_presence_absence.Rtab", sep="\t", index=False)
    (tmp_path / "gene_presence_absence.csv").write_text("")

    data = module.load_pangenome("panaroo", tmp_path)
    assert data.presence_absence.shape == (3, 3)
    assert data.family_metadata.shape[0] == 3
    assert data.family_metadata["annotation"].isna().all()
