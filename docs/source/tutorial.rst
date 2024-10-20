Tutorial
========

In this tutorial, we will obtain two binary complexes from the PDB, and model a Protac ternary complex (TC) with them. All files can be found in the tutorial folder in the `GitHub repo <https://example.com>`_.
We will replicate the first ever TC determined experimentally, which comprises the protein of interest (POI) BRD4bd2 bound to the Von Hippel Lindau (VHL) E3 ligase substrate receptor and bridged by the protac MZ1. First, we will run P4ward as if this structure has never been solved at all and make a model of it. Next, we will run P4ward in benchmarking mode, so that we can compare the results with the known crystal structure.

Obtaining the input files
------------------

First we need to obtain the files that P4ward needs to run.

- Obtain the binary complex of BRD4bd2 bound to MS417 inhibitor from the PDB, code `6DUV <https://www.rcsb.org/structure/6DUV>`_. We will name this file ``receptor_raw.pdb``
- Obtain the binary of VHL bound to ligand code 3JF, code `4W9H <https://www.rcsb.org/structure/4W9H>`_. This file will be ``ligase_raw.db``
- For each of these binary complexes, we also need to obtain a ``mol2`` file of their respective ligands. Therefore, we can scroll down on their PDB pages and download the instance coordinates for the ligands, choosing mol2 format.
    - For BRD4bd2, download the instance coordinates for the ligand mol2 `here <https://models.rcsb.org/v1/6duv/ligand?auth_seq_id=501&label_asym_id=C&encoding=mol2&filename=6duv_C_0S6.mol2>`_. We will save this file as ``receptor_ligand.mol2``
    - For VHL, find them `here <https://models.rcsb.org/v1/4w9h/ligand?auth_seq_id=301&label_asym_id=M&encoding=mol2&filename=4w9h_M_3JF.mol2>`_. We will save this as ``ligase_ligand.mol2``

Now, we need to clean up the protein files. Using the software of your choice, remove all except the protein's main chain of interest. For example, with ChimeraX we can:

- Open ``receptor_raw.pdb`` CONT HERE then save as ``receptor.pdb``
- Open ``ligase_raw.pdb`` CONT HERE then save as ``ligase.pdb``

Now, if we open ``receptor.pdb`` with ``receptor_ligand.mol2``, we should see the binary complex:
ADD IMG

And we should see the same for ``ligase.pdb`` and ``ligase_ligand.mol2``:
ADD IMG

Now we can obtain the Protacs. P4ward will take a smiles file with one or more smiles code per line. Each line is the code for an entire protac. For now, we can obtain the smiles code for protac MZ1 by copying it from its `PDB page <https://www.rcsb.org/ligand/759>`_ and we can paste it in a file called ``protacs.smiles``. Note that the smiles code can be followed by the molecule name. This way, the contents of ``protacs.smiles`` will be:

.. code-block:: text

   Cc1c(sc-2c1C(=N[C@H](c3n2c(nn3)C)CC(=O)NCCOCCOCCOCC(=O)N[C@H](C(=O)N4C[C@@H](C[C@H]4C(=O)NCc5ccc(cc5)C6=C(NCS6)C)O)C(C)(C)C)c7ccc(cc7)Cl)C mz1

An important step at this point is to examine the ligand files and the protac 2D structure determined in the smiles code and make sure they match, since P4ward will need to match them later. Bonds in aromatic rings, for example, can be represented in different ways, which could lead to incompatibilities. Visualize your smiles code in a program such as the web version of Marvin sketch, and the ligand files in a program that shows the bond orders, such as Pymol.
In this case, we see there is no difference between them:

ADD IMGS


Now we can add to our configuration file the names of the files we just prepared. Open a new file which we will call ``config.ini`` with the following contents:

.. code-block:: ini
   :caption: File: config.ini

   [general]
   receptor = receptor.pdb
   ligase = ligase.pdb
   receptor_ligand = receptor_ligand.mol2
   ligase_ligand = ligase_ligand.mol2


Checking the protac-ligand matches
------------------

P4ward also offers a simple way to check if the ligands and the protac match. If run the pipeline with the command:

.. code-block:: bash

   p4ward --config_file config.ini --check_lig_matches
