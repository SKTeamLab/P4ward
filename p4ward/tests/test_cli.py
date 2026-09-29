from pathlib import Path
import subprocess
import sys

# Repository root directory containing the p4ward package
REPO_ROOT = Path(__file__).resolve().parents[2]


def test_cli_help():
    """Verify that calling `python -m p4ward --help` works"""
    result = subprocess.run(
        [sys.executable, "-m", "p4ward", "--help"],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )

    # Check that the command ran successfully
    assert result.returncode == 0

    # Check that help text mentions usage and key CLI options
    assert "usage:" in result.stdout.lower()
    assert "--config_file" in result.stdout
    assert "--write_default" in result.stdout
    assert "--check_lig_matches" in result.stdout


def test_cli_short_help():
    """Verify that calling `python -m p4ward -h` works"""
    result = subprocess.run(
        [sys.executable, "-m", "p4ward", "-h"],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )

    assert result.returncode == 0
    assert "usage:" in result.stdout.lower()
