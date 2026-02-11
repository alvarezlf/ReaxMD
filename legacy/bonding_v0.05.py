#### Timer
from datetime import datetime
from timeit import default_timer as timer
t0=timer()                                          #timer
date0=datetime.now().strftime("%X %d/%m/%y")

import sys
import argparse
from math import sqrt
from itertools import islice
import os

#### PARAMETERS
cov_rad_file=os.path.expanduser("~/bin/covalent_radii.ssv") ### Path to the covalent radii database
r_alpha=1.3                                         ### Margin for comparing the bond distances vs r_cov

#### DATA                                   ####
parser=argparse.ArgumentParser()
parser.add_argument("trajectory_file", help="name of the trajectory file.")
parser.add_argument("a", help="A length of the cell.")
parser.add_argument("b", help="B length of the cell.")
parser.add_argument("c", help="C length of the cell.")
parser.add_argument("-o", "--options", type=str, help=("\"-o hmw\", where:\n"
                    "m writes the groups of molecules detected in the first step.\n"
                    "t writes a trajectory file \"trajectory_file.HSi.xyz\" with the atoms involved in reactions are replaced by Si.\n"
                    "w analyzes the Wannier centers."))
args=parser.parse_args()
filename=args.trajectory_file
A=float(args.a)
B=float(args.b)
C=float(args.c)
write_opt=args.options

fname=args.trajectory_file.replace(".xyz","")

with open(filename, 'r') as traj_file:
    traj_N = int(traj_file.readline().strip())
    traj_lines = sum(1 for i in traj_file)
traj_steps = int(traj_lines / (traj_N + 2))

with open(cov_rad_file, "r") as cov_rad_data:
    cov_rad = {line.split()[1]: float(line.split()[2]) for line in cov_rad_data if not line.startswith('#')}

#### DISTANCES FUNCTION
#"N\n
#comment\n
#A1 X1 Y1 Z1\n
#..."
def step2coor(file, N):
    file.readline()
    file.readline()
    coor = list(file.readline().split() for i in range(N))
    return (coor)
#[["A", "X", "Y", "Z"],...]
def traj2dist(coor_step, A, B, C):
    dist_step=list()
    for i in range(traj_N):
        dist_atom = [coor_step[i][0]]
        Xi, Yi, Zi = map(float, coor_step[i][1:4])
        for j in range(i, traj_N):
            Xj, Yj, Zj = map(float, coor_step[j][1:4])
            DXj=(abs(Xj-Xi)+A/2)%A-A/2
            DYj=(abs(Yj-Yi)+B/2)%B-B/2
            DZj=(abs(Zj-Zi)+C/2)%C-C/2
            Dij=sqrt(DXj**2 + DYj**2 + DZj**2)
            dist_atom.append(round(Dij, 3))
        dist_step.append(dist_atom)
    return (dist_step)
#[["A1", d11, d12, d13...],["A2",d22, d23...]...]
#### BONDING FUNCTION
#[["A", "X", "Y", "Z"],...]
def coor2bond(coor_step, A, B, C):
    bond_step=list()
    for i in range(traj_N):
        bond_atom=list([coor_step[i][0], i+1])
        Xi,Yi,Zi=map(float, coor_step[i][1:4])
        r_i=cov_rad[coor_step[i][0]]
        for j in range(i,traj_N):
            r_j=cov_rad[coor_step[j][0]]
            Xj,Yj,Zj=map(float, coor_step[j][1:4])
            DXj=(abs(Xj-Xi)+A/2)%A-A/2
            DYj=(abs(Yj-Yi)+B/2)%B-B/2
            DZj=(abs(Zj-Zi)+C/2)%C-C/2
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            bond_dist=(r_i+r_j)*r_alpha
            if (Dij <= bond_dist):
                bond_atom.append(j+1)
        bond_step.append(bond_atom)
    return(bond_step)
#[[A, ID, b1, b2...]...]
def bond2reac(b1, b2):
    if len(b1)!=len(b2):
        return("Bonding lists have different lengths")
    reac=list()
    for i in range(len(b1)):
        if b1[i]!=b2[i]:
            for bond in b1[i]:
                if (bond not in b2[i]):
                    AT1=b1[i][0]
                    ID1=b1[i][1]
                    AT2=b1[bond-1][0]
                    ID2=b1[bond-1][1]
                    reac.append(list([ID1, AT1, ID2, AT2, "broken"]))
            for bond in b2[i]:
                if (bond not in b1[i]):
                    AT1=b2[i][0]
                    ID1=b2[i][1]
                    AT2=b2[bond-1][0]
                    ID2=b2[bond-1][1]
                    reac.append(list([ID1, AT1, ID2, AT2, "formed"]))
    return(reac)
#[(Step, n1, a1, n2, a2, "formed/broken"),...]
t1=timer()                                          #timer

#### MAIN
if write_opt and ("m" in write_opt):
    print("Molecular structure analysis under construction.")
    exit()
if write_opt and ("w" in write_opt):
    print("Wannier centers analysis under construction.")
    exit()
    
with open(filename, 'r') as traj_file:
    coor_0=step2coor(traj_file,traj_N)
    bond_0=coor2bond(coor_0,A,B,C)
    rx_list=list()
    for s in range (1,traj_steps):
        coor_1=step2coor(traj_file,traj_N)
        bond_1=coor2bond(coor_1,A,B,C)
        rx=bond2reac(bond_0,bond_1)
        if rx:
            for i in rx:
                i.insert(0,s+1)
                print("Step {0}: Bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
            rx_list.append(i)
        coor_0=coor_1
        bond_0=bond_1
        sys.stdout.flush()
if (not rx_list):
    print("No reactions found!")
    exit()
t2=timer()                                          #timer
#### HIGHLIGHTED TRAJECTORY FILE
if (write_opt) and ("t" in write_opt) and (rx_list):
    atom_list=[[r[0], r[1], r[2]] for r in rx_list] + [[r[0], r[3], r[4]] for r in rx_list]
    atom_list.sort()
    with open(filename, "r") as traj_file, open(fname+"-HSi.xyz", "w") as htra_file:
        out_line=1
        atom_it=iter(atom_list)
        a = next(atom_it, None)
        for current_line, line in enumerate(traj_file, start=1):
            if a:
                atom_line=(a[0]-1)*(traj_N+2)+2+a[1]
                if (current_line == atom_line):
                    htra_file.write(line.replace(a[2], "Si"))
                    a=next(atom_it, None)
                else:
                    htra_file.write(line)
            else:
                htra_file.write(line)
            sys.stdout.flush()

t4=timer()                                          #timer
##########################################################
date1=datetime.now().strftime("%X %d/%m/%y")
def chrono(start,end):
    o_timer=end-start
    if o_timer>=60:
        return("{:.3f} minutes".format(o_timer/60))
    else:
        return("{:.3f} seconds".format(o_timer))

if (write_opt) and ("b" in write_opt):
    print("\n############################################################")
    print("##  STARTED ON       {0}".format(date0))
    print("##  FINISHED ON      {0}".format(date1))
    print("##  TRAJECTORY FILE  {0}".format(filename))
    print("##  COVALENT RADII   {0}".format(cov_rad_file))
    print("##  alpha            {0}".format(r_alpha))
    print("##  CONNECTIVITIES   ELAPSED TIME  {0}".format(chrono(t1,t2)))
    print("##  H. TRAJECTORY    ELAPSED TIME  {0}".format(chrono(t2,t4)))
    print("##  GLOBAL           ELAPSED TIME  {0}".format(chrono(t0,t4)))
    print("############################################################")






