Welcome to Protacs Pipeline's documentation!
============================================

Program Paths
-------------

* ``megadock``
* ``obabel``
* ``rxdock_root``


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

* ``cluster_poses``

Controls whether RMSD clustering should be performed for the protein poses. If ``filter_poses`` is activated, then the clustering is performed using only the poses that passed the filter. Using this flag alone does not make the clustering have any effect on the final ranking, make sure to also toggle ``rank_cluster_reps_only`` in the section ``protein_ranking``.
Default: ``False``

* ``clustering_cutoff``

RMSD cutoff for clustering, in angstroms.
Default: 2.0


Protein ranking
---------------

* ``final_ranking_megadock_score``

Use the original megadock score  as the final score for ranking the protein docked poses.
Default: True

* ``cluster_rep`` [``best`` | ``centroid``]

If clustering was performed, this flag chooses which protein pose, in each cluster, will represent the cluster. ``best`` means that the pose with the best score will be the cluster representative. ``centroid`` means the cluster's centroid will be its representative.
Default: ``best``

* ``rank_cluster_reps_only``

If clustering of the protein docking poses was performed, this flag controls if the final ranking should consider the clusters. If true, the ranking will include only cluster representatives (see flag ``cluster_rep``).
Default: False

* ``top_poses``

How many top protein poses will be forwarded to protac sampling.
Default: 10

* ``generate_poses`` [``all``|``top``| ``filtered``| ``filtered_centroids``|``top_centroids``]

The protein poses calculated by megadock are handled internally by the pipeline. Here, the user can choose which category of protein poses the pipeline should generate protein pdb files for.
Default: top

* ``generate_poses_altlocA``

Keep only alternate location A for when generating pdb files. Avoids some errors.
Default: True

* ``generated_poses_folder``

Name of the folder where the generated pdb files will reside.
Default: protein_docking


Linker sampling
---------------

* ``rdkit_sampling``

Use rdkit to perform protac sampling
Default: True

* ``protac_poses_folder``

Name of the folder where the generated protac sdf will reside.
Default: protac_sampling

* ``extend_flexible_small_linker``

If the linker consists of very few atoms, protac sampling will fail because small deviations on the extremities' positions will make bonds unfeasible. With this option, if the pipeline detects that the linker is short (see ``min_linker_length``), it will also consider more neighbouring atoms, from the ligands, as flexible (see ``extend neighbour number``).
Default: True

* ``extend_neighbour_number``

If ``extend_flexible_small_linker`` is turned on, then this flag controls how many neighbouring atoms should become flexible.
Default: 2

* ``min_linker_length``

If the protac's linker contains up to this many atoms, it is considered too short and can be extended if ``extend_flexible_small_linker`` is turned on.
Default: 2

* ``rdkit_number_of_confs``

How many protac poses to generate.
Default: 10

* ``rdkit_pose_rmsd_tolerance``

Some protac poses cannot be sampled while perfectly retaining the rigid ligands' positions. This flag controls how much deviation is allowed when this happens.
Default: 1.0 (angstroms)

* ``rdkit_time_tolerance``

Sometimes rdkit will get stuck for a very long time in a pose only to fail sampling. This flag sets a time limit to the time rdkit can spend in the sampling calculation for each pose. If the limit is reach, the pose is considered failed.
Default: 300 (seconds)

* ``extend_top_poses_sampled``

Extends how many protein poses are considered top (based on ``protein_ranking/top_poses``) so that the ``top_poses`` number of poses have successfully generated protac conformations. For example, if the user determined ``top_poses`` to be 10, then the top 10 protein poses will be forwarded to protac sampling. However, a few of these may not be optimal for protac conformation and so would fail at sampling. So the pipeline will try sampling for the 11th pose, 12th and so on, until exactly 10 poses have successfully generated protac conformations.
Default: True


Linker ranking
--------------

* ``protac_scoring_folder``

Name of the folder where the scored protac sdf files will reside.
Default: protac_scoring

* ``clash_detection``

Use biopython to detect if a protac conformation severely clashes with the proteins.
Default: True

* ``restrict_clash_to_linker``

Only consider the linker atoms, not the ligand ones, when looking for clashes.

* ``clash_threshold``

Distance below this value is considered a clash. Note: hydrogen atoms are not included.
Default: 2.5

* ``filter_clashed``

Turn on if the conformations with clashes should be considered unsuccessful.
Default: False

* ``max_clashes_allowed``

How many clashes are allowed before the pose is considered unsuccessful.
Default: 1

* ``rxdock_score``

Use RXdock for scoring the protac conformations.
Default: True

* ``rxdock_minimize``

Perform a a quick minimization with RXdock before scoring.
Default: True

* ``filter_scored_linkers``

Consider protac poses with positive scores unsuccessful.