Welcome to Protacs Pipeline's documentation!
============================================

Program Paths
-------------

* `megadock`
* `obabel`


General
-------

* ``overwrite`` = ``true`` | ``false``.

If there are records from previous runs of the pipeline, but they should be ignored and a overwritten with a new run starting from scratch. Default: ``false``.

* ``receptor``

Path to the receptor protein file, in pdb format. Default: ``receptor.pdb``.

* ``ligase``

Path to the ligase protein file, in pdb format. Default: ``ligase.pdb``.

* ``receptor_ligand``

Path to the receptor ligand file, in mol2 format. Default: ``receptor_ligand.mol2``.

* ``ligase_ligand``

Path to the ligase ligand file, in mol2 format. Default: ``ligase_ligand.mol2``.

* ``protac``

Path to the protac file, in smiles format. Default: ``protacs.smiles``


Megadock
--------

If Megadock is to be used for protein-protein docking (which is the only option right now) the following options regulate its usage:

* ``run_docking``

Toggles protein-protein docking with Megadock. Default: ``True``.

* ``num_predictions`` 

.. FILL THIS UP WITH DECENT INFORMATION
Default: ``2000``, from Megadock's default configuration values.

* ``num_predictions_per_rotation`` 

.. FILL THIS UP WITH DECENT INFORMATION
Default: ``1``, from Megadock's default configuration values.

* ``num_threads`` 

How many processor core to use for the megadock calculation.
Default: ``all``.

* ``num_threads`` 

How many processor core to use for the megadock calculation.
Default: ``all``.

* ``run_docking_output_file``

The main megadock output file path, where the docking results are placed. Default: ``megadock.out``

* ``run_docking_log_file``

The log file which will contain what megadock would normally print to the terminal. Default: ``megadock_run.log``

* ``filter_poses``

.. ADD REFERENCE TO EXPLANATION WHEN READY
Toggle to filter the protein poses based on the proximity of the ligands after docking. For more information on how this is done, please refer to [XXX]. The option can either be a number or ``auto`` to generate a proximity cutoff based on the protac size. 