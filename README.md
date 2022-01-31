# Analysis of protein nucleic-acid systems

This repository includes analysis scripts for RMSD, RMSF, distance and contact maps, and Watson-Crick hydrogen bond calculations for protein-nucleic acid molecular dynamics simulation trajectories. All scripts were written in Python and used MDAnalysis module.

*** Requirements:

Python versions 3+

MDAnalysis

*** Usage:

Python [options] filename.py

*** Examples:

python rmsd_mda.py --traj_path=trajectorydir --psf_path=psfdir --pdb_path=pdbdir --fit_selection='protein and (name CA or name CB)' --rmsd_selection='protein and name CA' --out_path=rmsd.dat

python rmsf_mda.py --traj_path=trajectorydir --psf_path=psfdir --pdb_path=pdbdir --fit_selection='protein and (name CA or name CB)' --rmsf_selection='protein and name CA' --out_path=rmsf.dat

python distance_mda.py --traj_path=trajectorydir --psf_path=psfdir --sela="protein" --selb="nucleic" --out=distance.dat

python contacts_mda.py --traj_path=trajectorydir --psf_path=psfdir --sela="protein" --selb="nucleic" --out=contact.dat

python hbonds_mda.py --traj_path=trajectorydir --psf_path=psfdir --out_path=hbond.dat

*** Citation

Adan Gallardo, Brandon Bogart, Bercem Dutagaci, Protein-Nucleic Acid Interactions for RNA Polymerase II Elongation Factors by Molecular Dynamics Simulations, BioRxiv (2022), doi: https://doi.org/10.1101/2022.01.28.478254

