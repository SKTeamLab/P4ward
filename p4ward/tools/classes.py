from pathlib import Path

class Protein:
    """
    Represent the main receptor target protein or the main ligase protein.
    """

    """
    attributes added/modified by the functions:
        md.fix_protein()
            - self.fixed_file
            - self.active_file
        md.minimize_proteins()
            - self.minim_file
            - self.active_file
        megadock.prep_structures()
            - self.mg_file
            - self.mg_file
        structure_tools.reduce()
            - self.reduced_file
    """

    
    def __init__(self, ptn_type, file, lig_file) -> None:
        """
        Starts off with the only 3 definitions required:
        pdb file with ligand bound, ligand chain ID, and ligand resnum
        automatically loads the biopython object for the protein
        """

        self.type = ptn_type # type is either 'receptor' or 'ligase'
        self.file = Path(file)
        self.lig_file = Path(lig_file)

        self.active_file = self.file # active file starts as the raw one

        if self.type == 'ligase':
            self.conformations = []
    

    def get_protein_struct(self, struct_attr='active_file'):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein. 
        """ 
        from ..run.structure_tools import load_biopython_structures
        structure_file = getattr(self, struct_attr)
        protein_struct = load_biopython_structures(structure_file=structure_file)
        return(protein_struct)


    def get_ligand_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein ligand. 
        """
        from ..run.structure_tools import load_biopython_structures

        ligand_struct = load_biopython_structures(
            structure_file=self.lig_file,
            mol2=True
        )
        return(ligand_struct)
    

    def active_confs(self):
        """
        Out of all the poses, return only the active ones
        """
        active_confs = [pose for pose in self.conformations if pose.active]
        return(active_confs)
    

    def get_triad_points(self):
        """
        Get the three points to represent the protein for clustering
        """
        import numpy as np

        a = self.get_ligand_struct().center_of_mass()
        b = self.get_protein_struct().center_of_mass()

        ab = b - a
        ab_normalized = ab / np.linalg.norm(ab)
        c = a + ab_normalized * 5  # 1 angstrom away from a

        dummy_vector = np.array([1, 0, 0])
        if (ab_normalized == dummy_vector).all():
            dummy_vector = np.array([0, 1, 0])  # Change dummy vector if it's parallel to ab

        perpendicular_vector = np.cross(ab_normalized, dummy_vector)
        perpendicular_vector_normalized = perpendicular_vector / np.linalg.norm(perpendicular_vector)
        d = c + perpendicular_vector_normalized * 5

        return(a, c, d)




class ProteinPose():
    """
    Represent a (docked) pose of the ligase protein.
    """

    """
    attributes added/modified by the functions:
        megadock.capture_scores()
            - self.active = True
            - self.megadock_score
        megadock.generate_poses()
            - self.file
        megadock.filter_poses()
            - self.active = Bool
            - self.filtered = Bool
        megadock.cluster()
            - self.rmsd
            - self.rmsd_reference
            - self.cluster
            - self.centroid = Bool
        megadock.zrank_rescore()
            - self.z_score
        rank.protein_poses()
            - self.active = Bool
            - self.top = Bool
        linker_sampling.rdkit_sampling()
            - self.protac_pos
        
    """

    def __init__(self, parent, pose_number) -> None:
        """
        requires its file and parent object. automatically generates
        the rest from the main ligase that is its "parent"
        """

        self.parent = parent
        self.pose_number = pose_number
        self.lig_file = parent.lig_file
        self.active = None
        self.file = None
        self.cluster_trend = None
        self.cluster_redund = None

        # add itself to the parent's conformations list
        if self not in parent.conformations:
            parent.conformations.append(self)
    

    def get_rotated_struct(self, struct_type, struct_attr='active_file'):
        """
        Use the rotate_atoms function to return a completely rotated
        coordinate set for the protein pose.
        """
        from Bio.PDB import Selection

        from ..run.megadock import rotate_atoms

        if struct_type == 'protein':
            ligase_obj = self.parent.get_protein_struct(struct_attr=struct_attr)
        elif struct_type == 'ligand':
            ligase_obj = self.parent.get_ligand_struct()
        else:
            raise Exception('Structure to capture must "protein" or "ligand".')
        ref_rotate = self.parent.rotate
        pose_rotate = self.rotate

        atoms = Selection.unfold_entities(ligase_obj, 'A')
        for atom in atoms:
            x,y,z, = atom.get_vector()
            newX, newY, newZ = rotate_atoms(
                (x, y, z),
                ref_rotation=ref_rotate,
                pose_rotation=pose_rotate
            )
            atom.set_coord((newX, newY, newZ))

        return(ligase_obj)


    def top_protac(self, protac_obj):
        protac_obj.protein_poses.append(self)




class Cluster():

    """
    attributes added/modified by the functions:
    """

    def __init__(self, clusterer, type) -> None:

        self.clusterer = clusterer
        # self.repr_centr = repr_centr
        # self.repr_best = repr_best
        self.type = type # 'redundancy', 'trend'
        self.clusters = {}
        # self.cutoff = self.clusterer.distance_threshold

    def get_all_confs(self):

        all_confs = []
        for i in self.clusters.values():
            all_confs.extend(i)

        return(all_confs)

    def get_centroid(self, cln):

        for pose_obj in self.clusters[cln]:
            if pose_obj in self.centroids:
                return(pose_obj)
    
    def get_cl_from_pose(self, pose_obj):

        for cln in self.clusters.keys():
            if pose_obj in self.clusters[cln]:
                return(cln)
        
    def get_cl_size(self, cln):

        return(len(self.clusters[cln]))
    



class Protac():

    """
    Represents a PLV.
    Each protac smiles code from the input file becomes a Protac obj.

    attributes added/modified by the functions:
        linker_sampling.rdkit_sampling()
            - self.indices_ligs
            - self.indices_link
        structure_tools.get_protac_dist_cuttoff()
            - self.dist_cutoff
        self.sample_unbound_confs()
            - self.unbound_confs
            - self.unbound_energy
    """

    def __init__(self, smiles, name, number) -> None:
        self.smiles = smiles
        self.name = name
        self.number = number
        self.poses = []
        self.protein_poses = []
        self.indices_ligs = None


    def active_poses(self):
        active_confs = [pose for pose in self.poses if pose.active]
        return(active_confs)
    
    def get_pose(self, pose_obj):
        for i in self.poses:
            if i.protein_parent is pose_obj:
                return(i)

    def sample_unbound_confs(self, num_unbound_confs=100):

        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import AllChem

        protac = Chem.MolFromSmiles(self.smiles)
        protac = Chem.AddHs(protac)
        params = AllChem.ETKDGv3()
        confs = AllChem.EmbedMultipleConfs(protac, numConfs=num_unbound_confs, params=params)
        ps = AllChem.MMFFGetMoleculeProperties(protac)

        # final energy value
        energies = []
        for confid in confs:
            ff = AllChem.MMFFGetMoleculeForceField(protac, ps, confId=confid)
            ff.Initialize()
            energies.append(ff.CalcEnergy())

        self.unbound_confs = protac
        self.num_confs = protac.GetNumConformers()
        self.unbound_energy = np.max(energies)
            

    
    def write_unbound_confs(self, num_unbound_confs=100, filename='unbound_protac.sdf'):

        from rdkit import Chem

        with open(filename, 'a+') as protac_file:
            for i in range(self.num_confs):
                molblock = Chem.MolToMolBlock(self.unbound_confs, confId=i, kekulize=False)
                protac_file.write(f'conf_{i}')
                protac_file.write(molblock)
                protac_file.write('$$$$\n')



class ProtacPose():

    """
    attributes added/modified by the functions:
        linker_sampling.rdkit_sampling()
            - self.active
            - self.file
        linker_scoring.rxdock_rescore()
            - self.scored_file
        rank.rescore()
            - self.rescore
    """

    def __init__(self, parent, protein_parent) -> None:
        self.parent = parent
        self.protein_parent = protein_parent
        self.linker_confs = []
        self.active = None

        self.protein_parent.protac_pose = self
        if self not in parent.poses:
            parent.poses.append(self)
    
    def active_confs(self):
        active_confs = [pose for pose in self.linker_confs if pose.active]
        return(active_confs)




class LinkerConf():

    """
    attributes added/modified by the functions:
        linker_sampling.rdkit_sampling()
        linker_sampling.capture_dock6_scores()
        linker_sampling.detect_clashes()
            - self.active
        linker_scoring.capture_rxdock_scores()
            - self.rx_score
        protac_scoring.get_unbound_rmsd()
            - self.unbound_stats
    """

    def __init__(self, parent, conf_number) -> None:
        self.parent = parent
        self.conf_number = conf_number

        if self not in parent.linker_confs:
            parent.linker_confs.append(self)
    