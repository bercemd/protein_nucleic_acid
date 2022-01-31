import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import numpy as np
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
    '--out_path',type=str, default='hbond.dat',
    help='Output path.')
arg = parser.parse_args()

def hbond_mda(u):
  #Watson Crick pair selection
  sela_a = "((resname ADE or resname DA) and name N1)"
  sela_t = "((resname THY or resname DT or resname DT3) and name O4)"
  sela_g = "((resname GUA or resname DG) and name O6)"
  sela_c = "((resname CYT or resname DC) and (name N3 or name O2))"

  selh_a = "((resname ADE or resname DA) and (name H61 or name H62))"
  selh_t = "((resname THY or resname DT or resname DT3) and name H3)"
  selh_g = "((resname GUA or resname DG) and (name H1 or name H22 or name H21))"
  selh_c = "((resname CYT or resname DC) and (name H41 or name H42))"

  hbonds_at = HBA(universe=u,acceptors_sel=sela_a,hydrogens_sel=selh_t)
  hbonds_ta = HBA(universe=u,acceptors_sel=sela_t,hydrogens_sel=selh_a)
  hbonds_gc = HBA(universe=u,acceptors_sel=sela_g,hydrogens_sel=selh_c)
  hbonds_cg = HBA(universe=u,acceptors_sel=sela_c,hydrogens_sel=selh_g)

  hbonds_at.run()
  count_at = hbonds_at.count_by_time()
  hbonds_ta.run()
  count_ta = hbonds_ta.count_by_time()
  hbonds_gc.run()
  count_gc = hbonds_gc.count_by_time()
  hbonds_cg.run()
  count_cg = hbonds_cg.count_by_time()

  return count_at, count_ta, count_gc, count_cg
dcds = []

dcds.append('%s/prod.1.dcd'%arg.traj_path)
dcds.append('%s/prod.2.dcd'%arg.traj_path)
dcds.append('%s/prod.3.dcd'%arg.traj_path)

U = mda.Universe(arg.psf_path, dcds, in_memory=True)
count_at, count_ta, count_gc, count_cg = hbond_mda(U)
outputfile = open(arg.out_path,"w")
for i in range(len(count_at)):
  print(i+1,count_at[i]+count_ta[i],count_gc[i]+count_cg[i],count_at[i]+count_ta[i]+count_gc[i]+count_cg[i],file=outputfile)
outputfile.close()

