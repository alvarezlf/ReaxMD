#### Timer
from datetime import datetime
from timeit import default_timer as timer
t0=timer()                                          #timer
date0=datetime.now().strftime("%X %d/%m/%y")

#### PARAMETERS
cov_rad_file="/home/alvarez/bin/covalent_radii.ssv" ### Absolute path to the covalent radii database
r_alpha=1.3                                         ### Margin for comparing the bond distances vs r_cov

#### DATA                                   ####
import argparse
from math import sqrt

parser=argparse.ArgumentParser()
parser.add_argument('filename')
parser.add_argument('a')
parser.add_argument('b')
parser.add_argument('c')
args=parser.parse_args()

A=float(args.a)
B=float(args.b)
C=float(args.c)

fname=args.filename.replace(".xyz","")

traj_file=open(args.filename, 'r')
traj_N=int(traj_file.readline().replace("\n", ""))
traj_lines=1
for i in traj_file:
    traj_lines+=1
traj_steps=traj_lines/int(traj_N+2)
traj_file.close()

cov_rad_data=open(cov_rad_file, "r")                #Covalent radii
cov_rad=list()
for l in cov_rad_data:
    line=l.split()
    if line[0]!="#":
        cov_rad.append(line)
cov_rad_data.close()

#### DISTANCES FUNCTION
#"N\n
#comment\n
#A1 X1 Y1 Z1\n
#..."
def step2coor(file, N):
    coor=list()
    file.readline()
    file.readline()
    for i in range(N):
        line=file.readline().split()
        coor.append(line)
    return (coor)
#[["A", "X", "Y", "Z"],...]
def traj2dist(coord_step, A, B, C):
    dist_step=list()
    for i in range(len(coord_step)):
        dist_atom=list()
        Xi=float(coord_step[i][1])
        Yi=float(coord_step[i][2])
        Zi=float(coord_step[i][3])
        dist_atom.append(coord_step[i][0])
        for j in range(i,len(coord_step)):
            DXj=abs(float(coord_step[j][1])-Xi)
            DYj=abs(float(coord_step[j][2])-Yi)
            DZj=abs(float(coord_step[j][3])-Zi)
            while DXj>A/2:
                DXj-=A
            while DYj>B/2:
                DYj-=B
            while DZj>C/2:
                DZj-=C
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            dist_atom.append(round(Dij,3))
        dist_step.append(dist_atom)
    return(dist_step)
#[["A1", d11, d12, d13...],["A2",d22, d23...]...]
#### BONDING FUNCTION
#[["A", "X", "Y", "Z"],...]
def coor2bond(coord_step, A, B, C):
    bond_step=list()
    for i in range(traj_N):
        bond_atom=list()
        bond_atom.append(coord_step[i][0])
        bond_atom.append(i+1)
        Xi=float(coord_step[i][1])
        Yi=float(coord_step[i][2])
        Zi=float(coord_step[i][3])
        for element in cov_rad:
            if (coord_step[i][0]==element[1]):
                r_i=float(element[2])
        for j in range(i,traj_N):
            for element in cov_rad:
                if (coord_step[j][0]==element[1]):
                    r_j=float(element[2])
            DXj=abs(float(coord_step[j][1])-Xi)
            DYj=abs(float(coord_step[j][2])-Yi)
            DZj=abs(float(coord_step[j][3])-Zi)
            while DXj>A/2:
                DXj-=A
            while DYj>B/2:
                DYj-=B
            while DZj>C/2:
                DZj-=C
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            bond_dist=(r_i+r_j)*r_alpha
            if (Dij<=bond_dist):
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
traj_file=open(args.filename, 'r')

#Step 1
coor_0=step2coor(traj_file,traj_N)
bond_0=coor2bond(coor_0,A,B,C)
#Steps 2 to the end
rx_list=list()
for s in range (1,traj_steps):
    coor_1=step2coor(traj_file,traj_N)
    bond_1=coor2bond(coor_1,A,B,C)
    rx=bond2reac(bond_0,bond_1)
    if rx!=[]:
        for i in rx:
            i.insert(0,s+1)
            print("Step {0}: Bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
        rx_list+=rx
    coor_0=coor_1
    bond_0=bond_1
if (len(rx_list)==0):
    print("No reactions found!\n")
    exit()
t2=timer()                                          #timer
    
#### HIGHLIGHTED TRAJECTORY FILE
if rx_list!=[]:
    atom_list=list()
    for r in rx_list:
        atom_list.append([r[0],r[1],r[2]])
        atom_list.append([r[0],r[3],r[4]])
    atom_list.sort()
#    print(rx_list)
#    print(atom_list)
    traj_file=open(args.filename, "r")
    htra_file_name=fname+"-HTRA.xyz"
    htra_file=open(htra_file_name, "w")
    out_step=1
    out_line=1
    for a in atom_list:
        atom_line=(a[0]-1)*(traj_N+2)+2+a[1]
        while (out_line!=atom_line):
            htra_file.write(traj_file.readline())
            out_line+=1
        htra_file.write(traj_file.readline().replace(a[2],"Si"))
        out_line+=1
    for i in traj_file:
        htra_file.write(i)
    traj_file.close()
    htra_file.close()

t4=timer()                                          #timer
##########################################################
date1=datetime.now().strftime("%X %d/%m/%y")
def chrono(start,end):
    o_timer=end-start
    if o_timer>=60:
        return("{:.3f} minutes".format(o_timer/60))
    else:
        return("{:.3f} seconds".format(o_timer))
               
print("\n############################################################")
print("##  STARTED ON       {0}".format(date0))
print("##  FINISHED ON      {0}".format(date1))
print("##  TRAJECTORY FILE  {0}".format(args.filename))
print("##  COVALENT RADII   {0}".format(cov_rad_file))
print("##  alpha            {0}".format(r_alpha))
print("##  CONNECTIVITIES   ELAPSED TIME  {0}".format(chrono(t1,t2)))
print("##  H. TRAJECTORY    ELAPSED TIME  {0}".format(chrono(t2,t4)))
print("##  GLOBAL           ELAPSED TIME  {0}".format(chrono(t0,t4)))
print("############################################################")






