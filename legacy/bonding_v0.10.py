import os
import sys
import argparse
from math import sqrt
from itertools import islice

#### INPUT                                   ####
parser=argparse.ArgumentParser()
parser.add_argument("trajectory_file", help="name of the trajectory file.")
parser.add_argument("a", type=float, help="A length of the cell in Angstrom.")
parser.add_argument("b", type=float, help="B length of the cell in Angstrom.")
parser.add_argument("c", type=float, help="C length of the cell in Angstrom.")
parser.add_argument("-a", "--alpha", type=float, help="Alpha factor for bond lengths as a*(r1+r2). Default value: 1.3")
parser.add_argument("-b", "--beta", type=float, help="Beta factor for Wannier centers proximities as b*(ri+g). Default value: 1.25")
parser.add_argument("-g", "--gamma", type=str, help="Gamma factor for Wannier centers proximities as b*(ri+g). Default value: 0.31")
parser.add_argument("-d", "--database", type=str, help="Alternative full path to covalent radii database. Default: /home/user/bin/covalent_radii.ssv")
parser.add_argument("-o", "--options", type=str, help=("\"-o mtw\", where:\n"
                    "m writes the groups of molecules detected in the first step.\n"
                    "t writes a trajectory file \"trajectory_file.HSi.xyz\" with the atoms involved in reactions are replaced by Si.\n"
                    "w analyzes the Wannier centers."))
parser.add_argument("-e", "--element", type=str, help=("\"-e El\" analyzes only the bonds involving element El.\n"
                    "\"-e E1-E2\" analyzes only bonds between elements E1 and E2."))
args=parser.parse_args()

#### INPUT
filename=args.trajectory_file
fname=args.trajectory_file.replace(".xyz","")
A=args.a
B=args.b
C=args.c
write_opt=args.options

#### Timer
benchmarking=(args.options) and ("b" in write_opt)      #timer
if benchmarking:
    from datetime import datetime
    from timeit import default_timer as timer
    t0=timer()
    date0=datetime.now().strftime("%X %d/%m/%y")
    def chrono(start,end):
        o_timer=end-start
        if o_timer>=60:
            return("{:.3f} minutes".format(o_timer/60))
        else:
            return("{:.3f} seconds".format(o_timer))

wc_activate=(args.options) and ("w" in write_opt)
th_activate=(args.options) and ("t" in write_opt)
if args.element:
    ELEM=args.element.split("-")
if args.database:
    cov_rad_file=args.database
else:
    cov_rad_file=os.path.expanduser("~/bin/covalent_radii.ssv")
if args.alpha:
    alpha=args.alpha
else:
    alpha=1.3
if args.beta:
    beta=args.beta
elif wc_activate and not args.beta:                 #try b,g=1.09,0.45 and b,g=1.095,0.39
    beta=1.25
if wc_activate and (not args.gamma):
    with open(cov_rad_file, "r") as cov_rad_data:
        #gamma={line.split()[1]: float(line.split()[3])*.01*2 for line in cov_rad_data if not line.startswith('#')}
        gamma={line.split()[1]: 0.31 for line in cov_rad_data if not line.startswith('#')}
if args.gamma:
    with open(cov_rad_file, "r") as cov_rad_data:
        gamma={line.split()[1]: float(args.gamma) for line in cov_rad_data if not line.startswith('#')}
if args.element and wc_activate:
    print ("Options \"-o w\" and \"-e\" are incompatible.")
with open(filename, 'r') as traj_file:
    traj_N = int(traj_file.readline().strip())
    traj_lines = sum(1 for i in traj_file)
traj_steps = int(traj_lines / (traj_N + 2))
with open(cov_rad_file, "r") as cov_rad_data:
    cov_rad = {line.split()[1]: float(line.split()[2]) for line in cov_rad_data if not line.startswith('#')}
#### DISTANCES FUNCTION
def step2coor(file, N):                     #file, 1 step -> coord. matrix
    file.readline()
    file.readline()
    coor=list()
    for i in range(1,N+1):
        line=file.readline().split()
        line.insert(0,i)
        coor.append(line)
    return (coor)
#[[ID1, "A1", "X1", "Y1", "Z1"],...]
def wannier(coor_w, coor_a,s):                #coord matrix -> atoms close to Wanniers
    wann_step=list()
    for i in coor_w:
        dist_atom = [i[:2]]
        Xi, Yi, Zi = map(float, i[2:5])
        for j in coor_a:
            r_j=cov_rad[j[1]]
            g_j=gamma[j[1]]
            Xj, Yj, Zj = map(float, j[2:5])
            DXj=(abs(Xj-Xi)+A/2)%A-A/2
            DYj=(abs(Yj-Yi)+B/2)%B-B/2
            DZj=(abs(Zj-Zi)+C/2)%C-C/2
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            if (Dij<=beta*(r_j+g_j)):
                dist_atom.append(j[:2])
        wann_step.append(dist_atom)
    return (wann_step)
#[[ID1, "A1", d11, d12, d13...],[ID2, "A2",d22, d23...]...]
#### BONDING FUNCTION
def coor2bond(coor_step, N):                #coordinates -> bonds
    bond_step=list()
    for i in range(N):
        bond_atom=list([coor_step[i][:2]])
        Xi,Yi,Zi=map(float, coor_step[i][2:5])
        r_i=cov_rad[coor_step[i][1]]
        for j in range(i,N):
            r_j=cov_rad[coor_step[j][1]]
            Xj,Yj,Zj=map(float, coor_step[j][2:5])
            DXj=(abs(Xj-Xi)+A/2)%A-A/2
            DYj=(abs(Yj-Yi)+B/2)%B-B/2
            DZj=(abs(Zj-Zi)+C/2)%C-C/2
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            bond_dist=(r_i+r_j)*alpha
            if (0<Dij<=bond_dist):
                bond_atom.append(coor_step[j][0:2])
        bond_step.append(bond_atom)
    return(bond_step)
#[[ID1, A1, [a1, id1], [a2, id2]...]...]
def bond2elem(c_elem, c_rest):              #coordinates -> bonds to specific element
    bond_step=list()
    for i in c_elem:
        bond_atom=list([i[:2]])
        Xi,Yi,Zi=map(float, i[2:5])
        r_i=cov_rad[i[1]]
        for j in c_rest:
            r_j=cov_rad[j[1]]
            Xj,Yj,Zj=map(float, j[2:5])
            DXj=(abs(Xj-Xi)+A/2)%A-A/2
            DYj=(abs(Yj-Yi)+B/2)%B-B/2
            DZj=(abs(Zj-Zi)+C/2)%C-C/2
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            bond_dist=(r_i+r_j)*alpha
            if (0<Dij<=bond_dist):
                bond_atom.append(j[0:2])
        bond_step.append(bond_atom)
    return(bond_step)
def bond2reac(b1, b2):                      #bonds -> reactions
    reac=list()
    for i in range(len(b1)):
        if b1[i]!=b2[i]:
            for bond in b1[i]:
                if (bond not in b2[i]):
                    AT1=b1[i][0][1]
                    ID1=b1[i][0][0]
                    AT2=bond[1]
                    ID2=bond[0]
                    reac.append(list([ID1, AT1, ID2, AT2, "broken"]))
            for bond in b2[i]:
                if (bond not in b1[i]):
                    AT1=b2[i][0][1]
                    ID1=b2[i][0][0]
                    AT2=bond[1]
                    ID2=bond[0]
                    reac.append(list([ID1, AT1, ID2, AT2, "formed"]))
    return(reac)
#[[Step, n1, a1, n2, a2, "formed/broken"],...]

def xyz_highlight(in_file, out_file, xx_list, Ael, Wel):
    atom_list1=[[r[0],r[1],r[2]] for r in xx_list]+[[r[0],r[3],r[4]] for r in xx_list]
    atom_list1.sort()
    atom_list=list()
    for i in atom_list1:
        if i not in atom_list:
            atom_list.append(i)
    with open(in_file, "r") as traj_file, open(out_file, "w") as htra_file:
        out_line=1
        atom_it=iter(atom_list)
        a=next(atom_it, None)
        for current_line, line in enumerate(traj_file, start=1):
            if a:
                atom_line=(a[0]-1)*(traj_N+2)+2+a[1]
                if (current_line == atom_line):
                    if a[2] != "X":
                        htra_file.write(line.replace(a[2], Ael))
                    elif a[2]=="X":
                        htra_file.write(line.replace(a[2], Wel))
                    a=next(atom_it, None)
                else:
                    htra_file.write(line)
            else:
                htra_file.write(line)
            sys.stdout.flush()

def checkdummy(coor):
    atoms=[i[1] for i in coor]
    if (not wc_activate) and ("X" in atoms):
        print("Error: Dummy atom \"X\" found in {0}!".format(filename))
        exit()
    if wc_activate and ("X" not in atoms):
        print("Error: No Wannier centers found in {0}!".format(filename))
        exit()
if benchmarking:                                        #timer
    t1=timer()

#### MAIN
#### El-Any bonds
#### E1-E2 bonds
rx_list=list()
wc_list=list()
if args.element:                            #-e
    with open(filename, 'r') as traj_file:
        coor_0=step2coor(traj_file,traj_N)
        checkdummy(coor_0)
        if len(ELEM)>2:
            print("Please specify a maximum of 2 elements")
            exit()
        elif len(ELEM)==2:
            coor_0=[i for i in coor_0 if i[1] in ELEM]
            bond_0=coor2bond(coor_0,len(coor_0))
        elif len(ELEM)==1:
            coor_0e=[i for i in coor_0 if i[1] in ELEM]
            coor_0r=[i for i in coor_0 if i[1] not in ELEM]
            bond_0=bond2elem(coor_0e,coor_0r)
        for s in range (1,traj_steps):
            coor_1=step2coor(traj_file,traj_N)
            if ELEM and len(ELEM)==2:
                coor_1=[i for i in coor_1 if i[1] in ELEM]
                bond_1=coor2bond(coor_1,len(coor_1))
            elif ELEM and len(ELEM)==1:
                coor_1e=[i for i in coor_1 if i[1] in ELEM]
                coor_1r=[i for i in coor_1 if i[1] not in ELEM]
                bond_1=bond2elem(coor_1e,coor_1r)
            rx=bond2reac(bond_0,bond_1)
            if rx:
                for i in rx:
                    i.insert(0,s+1)
                    print("Step {0}: Bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                rx_list.append(i)
            coor_0=coor_1
            bond_0=bond_1
            sys.stdout.flush()
if benchmarking:                                        #timer
    t2=timer()
#### Molecules
if args.options and ("m" in write_opt):
    print("Molecular structure analysis under construction.")
    exit()

#### Wannier centers
elif wc_activate:
    with open(filename, 'r') as traj_file:
        coor_0=step2coor(traj_file,traj_N)
        checkdummy(coor_0)
        coor_w0=[i for i in coor_0 if i[1]=="X"]
        coor_0=[i for i in coor_0 if i[1]!="X"]
        bond_0=coor2bond(coor_0,len(coor_0))
        wann_0=wannier(coor_w0, coor_0,0)
        for s in range (1,traj_steps):
            coor_1=step2coor(traj_file,traj_N)
            coor_w1=[i for i in coor_1 if i[1]=="X"]
            coor_1=[i for i in coor_1 if i[1]!="X"]
            bond_1=coor2bond(coor_1,len(coor_1))
            wann_1=wannier(coor_w0, coor_1,s)
            rx=bond2reac(bond_0,bond_1)
            if rx:
                for i in rx:
                    i.insert(0,s+1)
                    print("Step {0}: Bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                rx_list.append(i)
            wann=bond2reac(wann_0, wann_1)
            if wann:
                for i in wann:
                    i.insert(0,s+1)
                    if i[-1]=="formed":
                        i[-1]="approached"
                    if i[-1]=="broken":
                        i[-1]="separated from"
                    print("Step {0}: Wannier center {1} {5} {4}({3}).".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                    wc_list.append(i)
            coor_0=coor_1
            bond_0=bond_1
            coor_w0=coor_w1
            wann_0=wann_1
            sys.stdout.flush()
if benchmarking:                                        #timer
    t2=timer()
if (not args.element) and not wc_activate:
    with open(filename, 'r') as traj_file:
        coor_0=step2coor(traj_file, traj_N)
        checkdummy(coor_0)
        bond_0=coor2bond(coor_0, traj_N)
        for s in range (1,traj_steps):
            coor_1=step2coor(traj_file, traj_N)
            bond_1=coor2bond(coor_1, traj_N)
            rx=bond2reac(bond_0,bond_1)
            if rx:
                for i in rx:
                    i.insert(0,s+1)
                    print("Step {0}: Bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                rx_list.append(i)
            coor_0=coor_1
            bond_0=bond_1
            sys.stdout.flush()
if benchmarking:                                        #timer
    t2=timer()
if (not rx_list):
    print("No reactions found!")
if wc_activate and (not wc_list):
    print("No Wannier center's movement found!")

#### HIGHLIGHTED TRAJECTORY FILE
if th_activate:
    all_list=rx_list+wc_list
    if not all_list:
        print("Nothing to write!")
    else:
        xyz_highlight(filename, fname+"-HTr.xyz", all_list, "Si", "He")

if benchmarking:                                        #timer
    t3=timer()
    date1=datetime.now().strftime("%X %d/%m/%y")
    print("\n############################################################")
    print("##  STARTED ON       {0}".format(date0))
    print("##  FINISHED ON      {0}".format(date1))
    print("##  TRAJECTORY FILE  {0}".format(filename))
    print("##  COVALENT RADII   {0}".format(cov_rad_file))
    print("##  alpha            {0}".format(alpha))
    if wc_activate:
        print("##  beta             {0}".format(beta))
        print("##  gamma            {0}".format(gamma["H"]))
    print("##  REACTIONS        ELAPSED TIME  {0}".format(chrono(t1,t2)))
    if th_activate:
        print("##  H. TRAJECTORY    ELAPSED TIME  {0}".format(chrono(t2,t3)))
    print("##  GLOBAL           ELAPSED TIME  {0}".format(chrono(t0,t3)))
    print("############################################################")
