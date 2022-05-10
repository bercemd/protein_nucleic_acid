import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

import os,sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    '--traj_path',type=str, default='md.dcd',
    help='Trajectory path.')
parser.add_argument(
    '--psf_path',type=str, default='complex.psf',
    help='PSF path.')
parser.add_argument(
    '--out_path',type=str, default='base_distance.dat',
    help='Output path.')
arg = parser.parse_args()

def distance_mda(u,filename):
  select = "(segid DNAA)"
  sela = u.select_atoms(select, updating=True)
  selb = u.select_atoms(select, updating=True)
  sela_res = sela.residues
  selb_res = selb.residues

  i=0
  #iteration over the trajectory
  for frame in u.trajectory:
    print(frame)
    distance = ""
    for n in sela_res:
      for m in selb_res:
        if (n.resid-m.resid==1):
           xx = u.select_atoms("((segid DNAA) and (resid %s) and (name N1 or name C2 or name N3 or name C4 or name C5 or name C6 or name N7 or name C8 or name N9))"%n.resid)
           yy = u.select_atoms("((segid DNAA) and (resid %s) and (name N1 or name C2 or name N3 or name C4 or name C5 or name C6 or name N7 or name C8 or name N9))"%m.resid)
           xx_com = xx.center_of_mass()
           yy_com = yy.center_of_mass()
           distance = distance+str(np.linalg.norm(xx_com-yy_com))+" "
    print(distance,file=filename)
    i=i+1

dcds = []
dcds.append('%s/prod.1.dcd'%arg.traj_path)
dcds.append('%s/prod.2.dcd'%arg.traj_path)
dcds.append('%s/prod.3.dcd'%arg.traj_path)

U = mda.Universe(arg.psf_path, dcds, in_memory=True)
outputfile = open(arg.out_path,"w")
distance_mda(U,outputfile)
outputfile.close()

