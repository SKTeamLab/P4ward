

class Protein:
    """
    Represent the main receptor target protein or the main ligase protein.
    """

    """
    attributes added/modified by the functions:
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

        self.type = ptn_type #type is either 'receptor' or 'ligase'
        self.file = file
        self.lig_file = lig_file

        if self.type == 'ligase':
            self.conformations = []
    

    # def __getattr__(self, item):
    #     return(None)


    def get_protein_struct(self, struct_attr='file'):
        """
        Use the function load_biopython_structures to get the biopython object
        for the protein. 
        """ 
        from .structure_tools import load_biopython_structures
        structure_file = getattr(self, struct_attr)
        protein_struct = load_biopython_structures(structure_file=structure_file)
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


    # def __getattr__(self, item):
    #     return(None)
    

    def get_rotated_struct(self, struct_type, struct_attr='file'):
        """
        Use the rotate_atoms function to return a completely rotated
        coordinate set for the protein pose.
        """

        from ..run.megadock import rotate_atoms

        if struct_type == 'protein':
            ligase_obj = self.parent.get_protein_struct(struct_attr=struct_attr)
        elif struct_type == 'ligand':
            ligase_obj = self.parent.get_ligand_struct()
        else:
            raise Exception('Structure to capture must "protein" or "ligand".')
        ref_rotate = self.parent.rotate
        pose_rotate = self.rotate

        for model in ligase_obj:
            for chain in model:
                for res in chain:
                    for atom in res:
                        x,y,z, = atom.get_vector()
                        newX, newY, newZ = rotate_atoms(
                            (x, y, z),
                            ref_rotation=ref_rotate,
                            pose_rotation=pose_rotate
                        )
                        atom.set_coord((newX, newY, newZ))

        return(ligase_obj)