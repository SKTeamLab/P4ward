from ..tools import decorators

def get_interface(interface_struct, search_struct, cutoff=5.0):
    
    from Bio.PDB import NeighborSearch

    ns = NeighborSearch(list(search_struct.get_atoms()))
    interface = []

    for atom in interface_struct.get_atoms():
        if ns.search(atom.get_coord(), cutoff, level="R"):
            residue = atom.get_parent()
            if residue not in interface:
                interface.append(residue)
    
    return(interface)


def get_atom_pairs(interface_struct, search_struct, cutoff=5.0):

    from Bio.PDB import NeighborSearch

    ns = NeighborSearch(list(search_struct.get_atoms()))
    residue_pairs = []

    neighbours = []
    for residue in interface_struct.get_residues():
        neighbours = []
        for atom in residue.get_atoms():
            nss = ns.search(atom.get_coord(), cutoff, level="R")
            for ns_res in nss:
                if ns_res not in neighbours:
                    neighbours.append(ns_res)
        for neighbour in neighbours:
            pair = [residue, neighbour]
            if pair not in residue_pairs:
                residue_pairs.append(pair)

    return(residue_pairs)


def make_structure(*lists_of_residues, structure_name="structure"):

    from Bio.PDB import Structure, Model, Chain
    from string import ascii_uppercase as alc

    struct = Structure.Structure(structure_name)
    model = Model.Model(0)
    struct.add(model)

    for i in range(len(lists_of_residues)):

        residues = lists_of_residues[i]
        chain_id = alc[i]

        chain = Chain.Chain(chain_id)
        model.add(chain)
        for residue in residues:
            chain.add(residue)

    return(struct)



def calc_fnat(receptor_struct, ref_ligase_struct, pose_struct):

    ligase_if = get_atom_pairs(ref_ligase_struct, receptor_struct)
    rec_ligase_if = get_atom_pairs(receptor_struct, ref_ligase_struct)

    pose_if = get_atom_pairs(pose_struct, receptor_struct)
    rec_pose_if = get_atom_pairs(receptor_struct, pose_struct)

    ref_pairs = [*ligase_if, *rec_ligase_if]
    pose_pairs = [*pose_if, *rec_pose_if]
    correct_pairs = [i for i in pose_pairs if i in ref_pairs]

    fnat = len(correct_pairs) / len(ref_pairs)

    return(fnat)


def calc_lrms(receptor_struct, ref_ligase_struct, pose_struct):

    from copy import deepcopy
    from Bio.PDB import Superimposer
    from ..tools.structure_tools import get_rmsd

    if len(list(ref_ligase_struct.get_atoms())) <= len(list(receptor_struct.get_atoms())):
        l_rms = get_rmsd(ref_ligase_struct, pose_struct, backbone=True)
    
    else:
        rec_struct_align = deepcopy(receptor_struct)
        superimposer = Superimposer()
        superimposer.set_atoms(list(ref_ligase_struct.get_atoms()), list(pose_struct.get_atoms()))
        superimposer.apply(rec_struct_align)

        l_rms = get_rmsd(receptor_struct, rec_struct_align, backbone=True)
       
    return(l_rms)


def calc_irms(receptor_struct, ref_ligase_struct, ligase_obj, pose_obj):

    from copy import deepcopy
    from Bio.PDB import Selection
    from ..run.megadock import rotate_atoms
    from Bio.PDB import Superimposer
    from ..tools.structure_tools import get_rmsd

    # get residues in receptor that are within 10 of reflig
    rec_if = get_interface(
        interface_struct=receptor_struct,
        search_struct=ref_ligase_struct,
        cutoff=10.0
    )
    # get residues in reflig that are within 10 of receptor
    ref_lig_if = get_interface(
        search_struct=receptor_struct,
        interface_struct=ref_ligase_struct,
        cutoff=10.0
    )
    # build reference interface struct with both interfaces
    ref_if_struct = make_structure(rec_if, ref_lig_if, structure_name='ref_if_struct')

    # now we make a struct from interface residues of reflig
    lig_if_struct = make_structure(ref_lig_if)
    pose_if_struct = deepcopy(lig_if_struct)
    # and rotate it to match the pose_obj orientation
    atoms = Selection.unfold_entities(pose_if_struct, 'A')
    for atom in atoms:
        x,y,z, = atom.get_vector()
        newX, newY, newZ = rotate_atoms(
            (x, y, z),
            ref_rotation=ligase_obj.rotate,
            pose_rotation=pose_obj.rotate
        )
        atom.set_coord((newX, newY, newZ))

    pose_if = list(pose_if_struct.get_residues())
    # make the model interface struct with rec_if and pose_if
    model_if_struct = make_structure(rec_if, pose_if, structure_name='model_if_struct')

    # now we calculate the rmsd between them, but aligning coordinates first
    i_rms = get_rmsd(ref_if_struct, model_if_struct, fit=True, backbone=True)

    return(i_rms)


def calc_rank(fnat, lrms, irms):

    print(fnat, lrms, irms)

    if (
        ( fnat >= 0.5 ) and
        ( lrms <= 1.0 or irms <= 1.0 )
    ):
        rank = 'high'

    elif (
        ( fnat >= 0.3 and fnat < 0.5 ) and
        ( lrms <= 5.0 or irms <= 2.0 ) or
        ( fnat >= 0.5 and lrms > 1.0 and irms > 1.0 )
    ):
        rank = 'medium'

    elif (
        ( fnat >= 0.1 and fnat < 0.3 ) and
        ( lrms <= 10.0 or irms <= 4.0 ) or
        ( fnat >= 0.3 and lrms > 0.5 and irms > 2.0 )
    ):
        rank = 'acceptable'

    elif (
        ( fnat < 0.1 ) or
        ( lrms > 10.0 and irms > 4.0 )
    ):
        rank = 'incorrect'
    
    return(rank)


@decorators.user_choice
def benchmark(protac_objs, receptor_obj, ligase_obj, ref_ligase_file):

    from ..tools.structure_tools import load_biopython_structures

    for protac_obj in protac_objs:

        for pose_obj in protac_obj.protein_poses:

            receptor_struct = receptor_obj.get_protein_struct()
            ref_ligase_struct = load_biopython_structures(ref_ligase_file)
            pose_struct = pose_obj.get_rotated_struct('protein')

            fnat  = calc_fnat(receptor_struct, ref_ligase_struct, pose_struct)
            lrms = calc_lrms(receptor_struct, ref_ligase_struct, pose_struct)
            irms = calc_irms(receptor_struct, ref_ligase_struct, ligase_obj, pose_obj)
            rank = calc_rank(fnat, lrms, irms)

            pose_obj.capri = {'fnat':fnat, 'l_rms':lrms, 'i_rms':irms}
            pose_obj.capri_rank = rank
