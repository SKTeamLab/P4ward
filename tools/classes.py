

class Protein:
    """
    Represent the main receptor target protein or the main ligase protein.
    """

    """
    attributes added/modified by the functions:
        structure_tools.reduce()
            - self.reduced_file
        megadock.prep_structures()
            - self.prep_receptor_file
            - self.prep_ligase_file
    """

    
    def __init__(self, ptn_type, file, lig_file) -> None:
        """
        Starts off with the only 3 definitions required:
        pdb file with ligand bound, ligand chain ID, and ligand resnum
        automatically loads the biopython object for the protein
        """

        self.type = ptn_type #type is either 'receptor' or 'ligase'
        self.file = file
        self.lig_file = lig_file

        if self.type == 'ligase':
            self.conformations = []


    def get_protein_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein. 
        """ 
        from .structure_tools import load_biopython_structures
        protein_struct = load_biopython_structures(structure_file=self.file)
        return(protein_struct)


    def get_ligand_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein ligand. 
        """
        from .structure_tools import load_biopython_structures

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




class ProteinPose():
    """
    Represent a (docked) pose of the ligase protein.
    """

    """
    attributes added/modified by the functions:
        megadock.capture_scores()
            - self.megadock_score
            - self.active = True
        megadock.generate_poses()
            - self.file
        megadock.filter_poses()
            - self.active = Bool
            - self.file
        megadock.cluster()
            - self.rmsd
            - self.rmsd_reference
            - self.cluster
            - self.centroid = Bool
        megadock.zrank_rescore()
            - self.z_score
        linker_sampling.rdkit_sampling()
            - self.protac_file
        
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

        # add itself to the parent's conformations list
        if self not in parent.conformations:
            parent.conformations.append(self)


    def get_protein_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein. 
        """ 
        from .structure_tools import load_biopython_structures
        protein_struct = load_biopython_structures(structure_file=self.file)
        return(protein_struct)


    def get_ligand_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein ligand. 
        """
        from .structure_tools import load_biopython_structures

        ligand_struct = load_biopython_structures(
            structure_file=self.lig_file,
            mol2=True
        )
        return(ligand_struct)