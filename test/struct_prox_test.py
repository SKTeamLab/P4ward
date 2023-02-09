from protacs_pipeline.run import megadock
import os

CWD = '/home/paula/Documents/testing/protacs_pipeline/'
structures = os.listdir(os.path.join(CWD, 'protein_docking'))
structures = [os.path.join(CWD, 'protein_docking', i) for i in structures]

megadock.cluster(structures, 5.0, logger=None, choice=True)