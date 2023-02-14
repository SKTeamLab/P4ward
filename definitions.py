import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()
PIPELINE_NAME = 'protacs_pipeline'
PICKLE_FILE = os.path.join(CWD, PIPELINE_NAME+'-cpt.pickle')

