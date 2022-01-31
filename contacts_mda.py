import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.contacts import contact_matrix
from MDAnalysis.analysis.distances import distance_array

import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--traj_path',type=str, default='md',
    help='Trajectory path.')
parser.add_argument(
    '--psf_path',type=str, default='protein.psf',
    help='PSF path.')
parser.add_argument(
    '--sela',type=str, default='protein',
    help='First selection.')
parser.add_argument(
    '--selb',type=str, default='segid DNAA',
    help='Second selection.')
parser.add_argument(
    '--out_path',type=str, default='contacts.dat',
    help='Output path.')
arg = parser.parse_args()

def contacts_mda(u,selecta,selectb):
        #initial selection
        sel = u.select_atoms(selecta, updating=True)
        ref = u.select_atoms(selectb, updating=True)

        #getting residue names and residue numbers
        sel_res = sel.residues
        sel_resids = sel_res.resids
        ref_res = ref.residues
        ref_resids = ref_res.resids

        #empty array
        count_stack = np.zeros((len(sel_res), len(ref_res)))

        #iteration over the trajectory
        for frame in u.trajectory:
                print(frame)

                res_dist = np.full((len(sel_res), len(ref_res)), 100)

                #iterating through each of the contacting residues
                for n in sel_res:
                        for m in ref_res:
                                #selecting the atoms
                                sel_atoms = n.atoms.positions
                                ref_atoms = m.atoms.positions

                                #calculating the minimum distance between the residues
                                distance = distance_array(sel_atoms, ref_atoms, backend='OpenMP')
                                min_dist = np.min(distance)

                                #replacing the minimum distance
                                sel_i = (np.where(n == sel_res)[0][0])
                                ref_j = (np.where(m == ref_res)[0][0])
                                res_dist[sel_i, ref_j] = min_dist

                #cutoff distance
                cutoff=5
                contacts = contact_matrix(res_dist, cutoff)
                #Converting True to 1, and False to 0
                contact_int = contacts * 1
                #adding the contacts
                count_stack = count_stack + contact_int

        #Averaging the total number of contacts by the number of frames
        count_avg = count_stack / len(u.trajectory)
        counts = count_avg

        return counts, sel_resids, ref_resids

dcds = []
dcds.append('%s/prod.1.dcd'%arg.traj_path)
dcds.append('%s/prod.2.dcd'%arg.traj_path)
dcds.append('%s/prod.3.dcd'%arg.traj_path)

U = mda.Universe(arg.psf_path, dcds, in_memory=True)
contacts_arr, sela_resids, selb_resids = contacts_mda(U,arg.sela,arg.selb)
contacts_arr.shape
outputfile = open(arg.out_path,"w")
for i in range(len(contacts_arr)):
  for j in range(len(contacts_arr[i])):
    print(sela_resids[i],selb_resids[j],contacts_arr[i][j],file=outputfile)
outputfile.close()

