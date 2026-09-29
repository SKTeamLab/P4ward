from pathlib import Path
from p4ward.config import config
from p4ward.definitions import ROOT_DIR

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_make_config_defaults():
    """Verify that make_config loads default settings from default.ini"""
    # Load with tutorial config which specifies the files to use
    user_config_file = REPO_ROOT / "tutorial" / "config.ini"
    conf = config.make_config(user_config_file, ROOT_DIR)

    # Check that main sections are present
    assert conf.has_section("general")
    assert conf.has_section("megadock")
    assert conf.has_section("protein_filter")

    # Check that user values properly override defaults
    assert conf.get("general", "receptor") == "receptor.pdb"

    # Check that non-overridden defaults from default.ini are preserved
    assert conf.getint("general", "num_processors") > 0
