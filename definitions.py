import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()
PIPELINE_NAME = 'protacs_pipeline'
CPT_FILE = os.path.join(CWD, PIPELINE_NAME+'-cpt.pickle')
TRACKER_FILE = PIPELINE_NAME + '-tracker.pickle'

