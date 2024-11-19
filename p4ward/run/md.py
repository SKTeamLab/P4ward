import logging
from openmm.app import PDBFile
from ..tools import decorators

logger = logging.getLogger('p4ward')

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
                            hydrogens_only,
                            minimized_suffix='_minim'
):
    
    import openmm.app as omm
    import openmm.unit as omu
    from openmm import NoseHooverIntegrator, CustomExternalForce

    ff = omm.ForceField('amber14/protein.ff14SB.xml', 'implicit/gbn2.xml')
    
    for protein_obj in protein_objs:

        filepath = protein_obj.active_file
        minim_path = filepath.parent / (filepath.stem + minimized_suffix + filepath.suffix)

        # prepare system
        pdb = PDBFile(str(filepath))
        system = ff.createSystem(
            pdb.topology,
            nonbondedMethod=omm.NoCutoff
        )

        if hydrogens_only:

            force = 100000.0 * omu.kilocalories_per_mole/omu.angstroms**2
            restraint = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
            restraint.addGlobalParameter('k', force)
            restraint.addPerParticleParameter('x0')
            restraint.addPerParticleParameter('y0')
            restraint.addPerParticleParameter('z0')

            for atom_idx, atom in enumerate(pdb.topology.atoms()):
                if atom.element.symbol != 'H':
                    x0, y0, z0 = pdb.positions[atom_idx]
                    restraint.addParticle(atom_idx, [x0, y0, z0])
            
            system.addForce(restraint)

        # minimize
        integrator = NoseHooverIntegrator(300*omu.kelvin, 1/omu.picosecond, 0.002*omu.picoseconds)
        simulation = omm.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)

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


