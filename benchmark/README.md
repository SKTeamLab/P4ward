# Benchmarking 

This folder contains the structures that were used for benchmarking the pipeline as described [here](LINK). Users intending to extend the benchmarking studies or test different P4ward configurations with known TC structures can use these files to streamline their process.

## File description

There are 35 known TC structures included in this folder:

```text
5T35,6BN7,6BN8,6BN9,6BNB,6BOY,6HAX,6HAY,6HR2,6SIS,6ZHC,7JTO,7JTP,7KHH,7PI4,7Q2J,7S4E,7Z6L,7Z76,7Z77,7ZNT,8BB2,8BB4,8BB5,8BDS,8BDT,8BDX,8BEB,8FY0,8FY1,8FY2,8PC2,8QVU,8QW6,8QW7
```

In each of these, the files necessary to perform benchmarking with the bound proteins are:
- `receptor.pdb` - POI
- `ligase.pdb` - E3 substrate receptor. This file's rotation and translation were randomized to avoid biasing the protein docking calculations.
- `receptor_ligand.mol2`
- `ligase_ligand.mol2`
- `protac.smiles`
- `ref_ligase.pdb` - This is the binding pose of the E3 used for reference in benchmarking

The files for unbound protein benchmarking are:
- `receptor_unbound_opt.pdb` - POI
- `ligase_unbound_opt.pdb` - E3 substrate receptor. This file's rotation and translation were randomized to avoid biasing the protein docking calculations.
- `receptor_unbound_ligand.mol2`
- `ligase_unbound_ligand.mol2`
- `protac.smiles`
- `ref_ligase.pdb` - This is the binding pose of the E3 used for reference in benchmarking
