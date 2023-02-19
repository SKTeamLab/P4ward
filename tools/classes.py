

class Protein:
    """
    Represent the main receptor target protein or the main ligase protein.
    """

    """
    attributes added/modified by the functions:
    structure_tools.reduce()
        - self.reduced_file
    """

    
    def __init__(self, ptn_type, file, lig_chain, lig_resnum) -> None:
        """
        Starts off with the only 3 definitions required:
        pdb file with ligand bound, ligand chain ID, and ligand resnum
        automatically loads the biopython object for the protein
        """

        self.type = ptn_type #type is either 'receptor' or 'ligase'
        self.file = file
        self.lig_chain = lig_chain
        self.lig_resnum = lig_resnum

        if self.type == 'ligase':
            self.conformations = []


    def get_protein_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein. 
        """ 
        from .structure_tools import load_biopython_structures
        protein_struct = load_biopython_structures(protein=self.file)
        return(protein_struct)


    def get_ligand_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein ligand. 
        """
        from .structure_tools import load_biopython_structures

        _, ligand_struct = load_biopython_structures(
            protein=self.file,
            protein_ligand_chain=self.lig_chain,
            protein_ligand_resnum=self.lig_resnum
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
    """

    def __init__(self, parent, pose_number) -> None:
        """
        requires its file and parent object. automatically generates
        the rest from the main ligase that is its "parent"
        """

        self.parent = parent
        self.pose_number = pose_number
        self.lig_chain = parent.lig_chain
        self.lig_resnum = parent.lig_resnum
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
        protein_struct = load_biopython_structures(protein=self.file)
        return(protein_struct)


    def get_ligand_struct(self):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein ligand. 
        """
        from .structure_tools import load_biopython_structures

        _, ligand_struct = load_biopython_structures(
            protein=self.file,
            protein_ligand_chain=self.lig_chain,
            protein_ligand_resnum=self.lig_resnum
        )
        return(ligand_struct)
