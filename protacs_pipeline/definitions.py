from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent
CWD = Path.cwd()
PIPELINE_NAME = 'protacs_pipeline'
CPT_FILE = CWD / (PIPELINE_NAME+'-cpt.pickle')
TRACKER_FILE = CWD / (PIPELINE_NAME + '-tracker.pickle')

