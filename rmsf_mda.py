import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.rms as rms
from MDAnalysis.analysis.rms import RMSF
import MDAnalysis.analysis.align as align
from MDAnalysis.coordinates.memory import MemoryReader
import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--traj_path',type=str, default='md.dcd',
    help='Trajectory path.')
parser.add_argument(
    '--psf_path',type=str, default='protein.psf',
    help='PSF path.')
parser.add_argument(
    '--pdb_path',type=str, default='protein.pdb',
    help='PDB path.')
parser.add_argument(
    '--ref_selection',type=str, default='protein',
    help='Selection for reference.')
parser.add_argument(
    '--fit_selection',type=str, default='protein and (name CA or name CB)',
    help='Selection for alignment.')
parser.add_argument(
    '--rmsf_selection',type=str, default='protein and name CA',
    help='Selection for rmsf calculation.')
parser.add_argument(
    '--out_path',type=str, default='rmsf.dat',
    help='Output path.')
arg = parser.parse_args()


def rmsf_mda(u_in_mem, pdb, selection=arg.ref_selection):
	protein = u_in_mem.select_atoms(selection)

	#fitting trajectory to initial frame
	init_frame = mda.Universe(pdb)
	prealigner = align.AlignTraj(u_in_mem, init_frame, select=arg.fit_selection, in_memory=True).run()

	#building average structure
	reference_coordinates = u_in_mem.trajectory.timeseries(asel=protein).mean(axis=1)
	reference = mda.Merge(protein).load_new(reference_coordinates[:, None, :], order="afc")

	#align whole trajectory to average structure
	aligner = align.AlignTraj(u_in_mem, reference, select=arg.fit_selection, in_memory=True).run()

	mobile =  u_in_mem.select_atoms(arg.rmsf_selection)
	resnums = mobile.resnums
	rmsfer = RMSF(mobile).run(start=0, stop=(len(u_in_mem.trajectory)+1), step=1)
	return rmsfer.rmsf, resnums

U = mda.Universe(arg.psf_path, arg.traj_path, in_memory=True)	

outputfile = open(arg.out_path,"w")
data,  r_nums = rmsf_mda(U, arg.pdb_path, selection=arg.ref_selection)
for i in range(len(data)):
  print(data[i], r_nums[i],file=outputfile)
outputfile.close()
