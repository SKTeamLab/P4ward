from openmm.app import PDBFile, ForceField
from pdbfixer.pdbfixer import PDBFixer
from ..tools import decorators
from ..tools.logger import logger


@decorators.user_choice
@decorators.track_run
def fix_proteins(*protein_objs, fixed_suffix='_fixed', ignore_extremities=True, ph=7.0):
    """
    accepts any number of protein_obj (Path obj) to fix using pdbfixer
    creates attribute with fixed filepath and makes it active
    """

    for protein_obj in protein_objs:

        filepath = protein_obj.active_file
        fixed_path = filepath.parent / (filepath.stem + fixed_suffix + filepath.suffix)

        fixer = PDBFixer(str(filepath))

        fixer.findMissingResidues()
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()

        if ignore_extremities:
            delkeys = []
            for key in keys:
                chain = chains[key[0]]
                if key[1] == 0 or key[1] == len(list(chain.residues())):
                    delkeys.append(key)
            for key in delkeys:
                del(fixer.missingResidues[key])

        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(ph)

        with open(fixed_path, 'w+') as fixedfile:
            PDBFile.writeFile(fixer.topology, fixer.positions, fixedfile)

        # add fixed file to attributes and make it the active file
        protein_obj.fixed_file = fixed_path
        protein_obj.active_file = protein_obj.fixed_file

    # logit
    logger.info('Fixed proteins and saved them as: '+ ', '.join([str(i.fixed_file) for i in protein_objs]))


@decorators.user_choice
@decorators.track_run
def get_protein_charges(*protein_objs):

    from openmm import NonbondedForce
    ff = ForceField('amber14/protein.ff14SB.xml')

    for protein_obj in protein_objs:

        pdb = PDBFile(str(protein_obj.active_file)) # ideally should be fixed_file
        system = ff.createSystem(pdb.topology)

        nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
        charges = []
        for i in range(system.getNumParticles()):
            charge, _, _ = nonbonded.getParticleParameters(i)
            charges.append(charge._value)
        
        protein_obj.charges = charges


