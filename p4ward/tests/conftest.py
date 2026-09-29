import sys
from pathlib import Path

# Ensure the repository root is on sys.path so 'p4ward' can always be imported
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
