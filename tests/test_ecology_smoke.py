from pathlib import Path
import importlib.util
import sys


def load_module():
    scripts = Path(__file__).parents[1] / "scripts"
    sys.path.insert(0, str(scripts))
    path = scripts / "pangenome_ecology.py"
    spec = importlib.util.spec_from_file_location("pangenome_ecology", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_pairwise_metrics():
    module = load_module()
    result = module.pairwise_set_metrics(
        {"a": {"x", "y"}, "b": {"y", "z"}},
        {"a": 2, "b": 2},
    )
    assert len(result) == 1
    assert result.loc[0, "shared_families"] == 1
    assert result.loc[0, "jaccard"] == 1 / 3
