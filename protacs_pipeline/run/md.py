from openmm.app import PDBFile
from ..tools import decorators
from ..tools.logger import logger


@decorators.user_choice
@decorators.track_run
def fix_proteins(*protein_objs, fixed_suffix='_fixed', ignore_extremities=True, ph=7.0):
    """
    accepts any number of protein_obj (Path obj) to fix using pdbfixer
    creates attribute with fixed filepath and makes it active
    """
    from pdbfixer.pdbfixer import PDBFixer

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
def minimize_proteins(
                            *protein_objs,
                            maxiterations=0,
                            minimized_suffix='_minim'
):
    
    import openmm.app as omm
    import openmm.unit as omu
    from openmm import NoseHooverIntegrator

    ff = omm.ForceField('amber14/protein.ff14SB.xml', 'implicit/gbn2.xml')
    
    for protein_obj in protein_objs:

        filepath = protein_obj.active_file
        minim_path = filepath.parent / (filepath.stem + minimized_suffix + filepath.suffix)

        # prepare system
        pdb = PDBFile(str(filepath))
        system = ff.createSystem(
            pdb.topology,
            nonbondedMethod=omm.NoCutoff,
            constraints=omm.HBonds
        )
        integrator = NoseHooverIntegrator(300*omu.kelvin, 1/omu.picosecond, 0.002*omu.picoseconds)
        simulation = omm.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)

        # minimize
        logger.info(f'Minimizing energy for {filepath.name} ...')
        simulation.minimizeEnergy(maxIterations=maxiterations)

        # save pdb file
        minimized_positions = simulation.context.getState(getPositions=True).getPositions()
        with open(minim_path, 'w+') as minimfile:
            omm.PDBFile.writeFile(simulation.topology, minimized_positions, file=minimfile)
        
        # add the attribute file path and activate it
        protein_obj.minim_file = minim_path
        protein_obj.active_file = protein_obj.minim_file

        logger.info(f'Saved minimized file as {protein_obj.minim_file}')



# @decorators.user_choice
# @decorators.track_run
def get_protein_charges(*protein_objs):

    from openmm import NonbondedForce
    from openmm.app import ForceField

    ff = ForceField('amber14/protein.ff14SB.xml')

    for protein_obj in protein_objs:

        pdb = PDBFile(str(protein_obj.active_file)) # ideally should be fixed
        system = ff.createSystem(pdb.topology)

        nonbonded = [f for f in system.getForces() if isinstance(f, NonbondedForce)][0]
        charges = []
        for i in range(system.getNumParticles()):
            charge, _, _ = nonbonded.getParticleParameters(i)
            coords = ( # mult by 10 because openmm uses nm
                round(pdb.positions[i]._value.x * 10, 5),
                round(pdb.positions[i]._value.y * 10, 5),
                round(pdb.positions[i]._value.z * 10, 5),
            )
            charges.append((coords, charge._value))
        
        protein_obj.charges = charges


