#### bonding.py
#### Writen by L. Alvarez~
#### Covalent radii by B. Cordero et. al., Dalton Transactions, 2008, 21. https://doi.org/10.1039/B801115J
version=0.13

import os
import sys
import argparse
from math import sqrt
from itertools import islice
from collections import Counter

#### INPUT
parser=argparse.ArgumentParser()
parser.add_argument("trajectory_file", help="name of the trajectory file.")
parser.add_argument("a", type=float, help="A length of the cell in Angstrom.")
parser.add_argument("b", type=float, help="B length of the cell in Angstrom.")
parser.add_argument("c", type=float, help="C length of the cell in Angstrom.")
parser.add_argument("-a", "--alpha", type=float, help="Alpha value for bond lengths as a*(r1+r2). Default value: 1.3")
parser.add_argument("-b", "--beta", type=float, help="Beta value for Wannier centers proximities as b*(ri+g). Default value: 1.25")
parser.add_argument("-g", "--gamma", type=float, help="Gamma value for Wannier centers proximities as b*(ri+g). Default value: 0.31")
parser.add_argument("-r", "--radii", type=str, help="Alternative full path to covalent radii data. Default: /home/user/bin/bonding/covalent_radii.ssv")
parser.add_argument("-o", "--options", type=str, help=("\"-o mtwf\", where:\n"
                    "m writes a file \"trajectory-molecules.txt\" with the reactive molecules.\n"
                    "t writes a trajectory file \"trajectory_file-H.xyz\", where the atoms involved in reactions are replaced by Si and the Wannier centers (if wt or tw ) by He.\n"
                    "w analyzes the Wannier centers.\n"
                    "f filters the Wannier centers close to reactive atoms. Use with option \"w\""))
parser.add_argument("-e", "--element", type=str, help=("\"-e El\" analyzes only the bonds involving element El.\n"
                    "\"-e E1-E2\" analyzes only bonds between elements E1 and E2."))
args=parser.parse_args()

filename=args.trajectory_file
fname=args.trajectory_file.replace(".xyz","")
A=args.a
B=args.b
C=args.c
write_opt=args.options

#### CONFIGURATION
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
wc_filter=wc_activate and ("f" in write_opt)
th_activate=(args.options) and ("t" in write_opt)
if args.element:
    ELEM=args.element.split("-")
if args.radii:
    with open(args.radii, "r") as cov_rad_data:
        cov_rad = {line.split()[1]: float(line.split()[2]) for line in cov_rad_data if not line.startswith('#')}
else:
    cov_rad={
    'H': 0.31,  'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,  'C': 0.76,  
    'N': 0.71,  'O': 0.66,  'F': 0.57,  'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 
    'Al': 1.21, 'Si': 1.11, 'P': 1.07,  'S': 1.05,  'Cl': 1.02, 'Ar': 1.06, 
    'K': 2.03,  'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53,  'Cr': 1.39, 
    'Mn': 1.61, 'Fe': 1.52, 'Co': 1.50, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 
    'Ga': 1.22, 'Ge': 1.20, 'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16, 
    'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90,  'Zr': 1.75, 'Nb': 1.64, 'Mo': 1.54, 
    'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44, 
    'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39,  'Xe': 1.40, 
    'Cs': 2.44, 'Ba': 2.15, 'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 
    'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96, 'Tb': 1.94, 'Dy': 1.92, 
    'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75, 
    'Ta': 1.70, 'W': 1.62,  'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 
    'Au': 1.36, 'Hg': 1.32, 'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 
    'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21, 'Ac': 2.15, 'Th': 2.06, 
    'Pa': 2.00, 'U': 1.96,  'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69
    }
if args.alpha:
    alpha=args.alpha
else:
    alpha=1.3
if args.beta:
    beta=args.beta
elif wc_activate and not args.beta:
    beta=1.25
if wc_activate and (not args.gamma):
    gamma=0.31
if args.gamma:
    gamma=args.gamma
if args.element and wc_activate:
    print ("Options \"-o w\" and \"-e\" are incompatible.")
with open(filename, 'r') as traj_file:
    traj_N = int(traj_file.readline().strip())
    traj_lines = sum(1 for i in traj_file)+1
traj_steps = int(traj_lines/(traj_N + 2))

#### FUNCTIONS
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
def wannier(coor_w, coor_a,s):              #coord matrix -> atoms close to Wanniers
    wann_step=list()
    for i in coor_w:
        dist_atom = [i[:2]]
        Xi, Yi, Zi = map(float, i[2:5])
        for j in coor_a:
            r_j=cov_rad[j[1]]
            Xj, Yj, Zj = map(float, j[2:5])
            DXj=(abs(Xj-Xi)+A/2)%A-A/2
            DYj=(abs(Yj-Yi)+B/2)%B-B/2
            DZj=(abs(Zj-Zi)+C/2)%C-C/2
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            if (Dij<=beta*(r_j+gamma)):
                dist_atom.append(j[:2])
        wann_step.append(dist_atom)
    return (wann_step)
#[[ID1, "A1", d11, d12, d13...],[ID2, "A2",d22, d23...]...]
#### BONDING FUNCTION
def coor2bond(coor_step, N):                
    bond_step = []
    for i in range(N):
        Xi,Yi,Zi=map(float, coor_step[i][2:5])
        r_i=cov_rad[coor_step[i][1]]
        bond_atom=[tuple(coor_step[i][:2])]
        for j in range(i+1, N):
            r_j=cov_rad[coor_step[j][1]]
            Xj,Yj,Zj=map(float, coor_step[j][2:5])
            DXj=(Xj-Xi+A/2)%A-A/2
            DYj=(Yj-Yi+B/2)%B-B/2
            DZj=(Zj-Zi+C/2)%C-C/2
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            bond_dist=(r_i+r_j)*alpha
            if Dij <= bond_dist:
                bond_atom.append(tuple(coor_step[j][:2]))
        bond_step.append(bond_atom)
    return (bond_step)
#[[ID1, A1, [a1, id1], [a2, id2]...]...]
def bond2elem(c_elem, c_rest):
    bond_step = []
    for i in c_elem:
        Xi,Yi,Zi=map(float, i[2:5])
        r_i=cov_rad[i[1]]
        bond_atom=[tuple(i[:2])]
        for j in c_rest:
            r_j=cov_rad[j[1]]
            Xj,Yj,Zj=map(float, j[2:5])
            DXj=(Xj-Xi+A/2)%A-A/2
            DYj=(Yj-Yi+B/2)%B-B/2
            DZj=(Zj-Zi+C/2)%C-C/2
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            bond_dist=(r_i+r_j)*alpha
            if Dij <= bond_dist:
                bond_atom.append(tuple(j[:2]))
        bond_step.append(bond_atom)
    return (bond_step)
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
def join_molecules(bond_list):
    molecules=[]
    for g in bond_list:
        g=set(g)
        for m in molecules:
            if not g.isdisjoint(m):
                m.update(g)
                break
        else:
            molecules.append(g)
    check=True
    while check:
        check=False
        for i in range(len(molecules)):
            for j in range(i + 1, len(molecules)):
                if not molecules[i].isdisjoint(molecules[j]):
                    molecules[i].update(molecules[j])
                    del molecules[j]
                    check = True
                    break
            if check:
                break
    mol_out=[sorted(i) for i in molecules]
    return (mol_out)
def mol_formula(mol_list):
    form_dict = {}
    for m in mol_list:
        atom_count = Counter(a[1] for a in m)
        formula = ''
        if 'C' in atom_count:
            formula += 'C{}'.format(atom_count['C'] if atom_count['C'] > 1 else '')
            if 'H' in atom_count:
                formula += 'H{}'.format(atom_count['H'] if atom_count['H'] > 1 else '')
            for atom in sorted(atom_count):
                if atom != 'C' and atom != 'H':
                    count = atom_count[atom]
                    formula += '{}{}'.format(atom, count if count > 1 else '')
        else:
            for atom in sorted(atom_count):
                count = atom_count[atom]
                formula += '{}{}'.format(atom, count if count > 1 else '')
        if formula in form_dict:
            form_dict[formula] += 1
        else:
            form_dict[formula] = 1
    return (form_dict)

def reactive_molecules(list_r,list_p):
    list_rx=[i for i in list_r if i not in list_p]
    return(list_rx)

if benchmarking:                                        #
    t1=timer()

#### MAIN
rx_list=list()
wc_list=list()
if args.element:                                        #-e
    with open(filename, 'r') as traj_file:
        coor_0=step2coor(traj_file,traj_N)
        checkdummy(coor_0)
        if len(ELEM)>2 or not ELEM:
            print("Please specify 1 or 2 elements as \"-e El\" or \"-e E1-E2\"")
            exit()
        elif len(ELEM)==2:
            coor_0a1=[i for i in coor_0 if i[1] == ELEM[0]]
            coor_0a2=[i for i in coor_0 if i[1] == ELEM[1]]
            bond_0=bond2elem(coor_0a1,coor_0a2)
            bond_0_save=bond_0
        elif len(ELEM)==1:
            coor_0e=[i for i in coor_0 if i[1] in ELEM]
            coor_0r=[i for i in coor_0 if i[1] not in ELEM]
            bond_0=bond2elem(coor_0e,coor_0r)
            bond_0_save=bond_0
        for s in range (1,traj_steps):
            coor_1=step2coor(traj_file,traj_N)
            if ELEM and len(ELEM)==2:
                coor_1a1=[i for i in coor_1 if i[1]==ELEM[0]]
                coor_1a2=[i for i in coor_1 if i[1]==ELEM[1]]
                bond_1=bond2elem(coor_1a1,coor_1a2)
                rx_1=bond2reac(bond_0,bond_1)
                rx=[]
                for i in rx_1:
                    if [i[2],i[3],i[0],i[1],i[4]] not in rx:
                        rx.append(i)
            elif ELEM and len(ELEM)==1:
                coor_1e=[i for i in coor_1 if i[1] in ELEM]
                coor_1r=[i for i in coor_1 if i[1] not in ELEM]
                bond_1=bond2elem(coor_1e,coor_1r)
                rx=bond2reac(bond_0,bond_1)
            if rx:
                for i in rx:
                    i.insert(0,s+1)
                    print("step {0}: bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                rx_list.append(i)
            coor_0=coor_1
            bond_0=bond_1
            sys.stdout.flush()
    bond_S_save=bond_1
elif wc_activate:                                       #-o w
    with open(filename, 'r') as traj_file:
        coor_0=step2coor(traj_file,traj_N)
        checkdummy(coor_0)
        coor_w0=[i for i in coor_0 if i[1]=="X"]
        coor_0=[i for i in coor_0 if i[1]!="X"]
        bond_0=coor2bond(coor_0,len(coor_0))
        bond_0_save=bond_0
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
                    print("step {0}: bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                    rx_list.append(i)
            wann=bond2reac(wann_0, wann_1)
            if wann:
                for i in wann:
                    i.insert(0,s+1)
                    if i[-1]=="formed":
                        i[-1]="approached"
                    if i[-1]=="broken":
                        i[-1]="separated from"
                    #if not wc_filter:
                        #print("Step {0}: Wannier center {1} {5} {4}({3}).".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                    wc_list.append(i)
            coor_0=coor_1
            bond_0=bond_1
            coor_w0=coor_w1
            wann_0=wann_1
            sys.stdout.flush()
    bond_S_save=bond_1

if not wc_filter:
    for i in wc_list:
        print("step {0}: Wannier center {1} {5} {4}({3}).".format(i[0], i[1], i[2], i[3], i[4], i[5]))
if wc_filter:                                           #-o f
    rx_atom=set([i[1] for i in rx_list]+[i[3] for i in rx_list])
    wc_f=[w for w in wc_list if w[3] in rx_atom]
    for i in wc_f:
        print("step {0}: Wannier center {1} {5} {4}({3}).".format(i[0], i[1], i[2], i[3], i[4], i[5]))
if (not args.element) and not wc_activate:
    with open(filename, 'r') as traj_file:
        coor_0=step2coor(traj_file, traj_N)
        checkdummy(coor_0)
        bond_0=coor2bond(coor_0, traj_N)
        bond_0_save=bond_0
        for s in range (1,traj_steps):
            coor_1=step2coor(traj_file, traj_N)
            bond_1=coor2bond(coor_1, traj_N)
            rx=bond2reac(bond_0,bond_1)
            if rx:
                for i in rx:
                    i.insert(0,s+1)
                    print("step {0}: bond ({1}){2}-{4}({3}) is {5}.".format(i[0], i[1], i[2], i[3], i[4], i[5]))
                rx_list.append(i)
            coor_0=coor_1
            bond_0=bond_1
            sys.stdout.flush()
    bond_S_save=bond_1

#### MOLECULES
#if args.options and ("m" in write_opt):
    #print("Molecular structure analysis under construction.")
    #exit()
if benchmarking:                                        #
    t2=timer()

#### REACTION
if rx_list:
    mole_0=join_molecules(bond_0_save)
    form_0=mol_formula(mole_0)
    mole_S=join_molecules(bond_S_save)
    form_S=mol_formula(mole_S)
    for i in form_0:
        if i not in form_S:
            form_S[i]=0
    for i in form_S:
        if i not in form_0:
            form_0[i]=0
    form_Rx={}
    Reac=''
    Prod=''
    for i in form_S:
        form_Rx[i]=form_S[i]-form_0[i]
    for i in form_Rx:
        if form_Rx[i]>0:
            Prod+=' + {0} {1}'.format(form_Rx[i], i)
            Prod=Prod.replace('1 ','')
        elif form_Rx[i]<0:
            Reac+='{0} {1} + '.format(-form_Rx[i], i)
            Reac=Reac.replace('1 ','')
    net_rx=[i for i in form_Rx if form_Rx[i]!=0]
    if net_rx:
        print('\nNET REACTION:')
        print(Reac[:-2]+' -> '+Prod[2:])
    else:
        print('\nNO NET REACTION FOUND')
    
if args.options and ('m' in write_opt) and rx_list:     #-o m
    react_list=reactive_molecules(mole_0, mole_S)
    with open(fname+'-molecules.txt', 'w') as mole_file:
        mole_file.write('{0}\n{1}\n'.format(filename,traj_N))
        mole_file.write(r'# List of reactive atoms generated automatically by bonding.py'+'\n')
        for i in react_list:
            for j in i:
                mole_file.write(str(j[0])+'\n')
            mole_file.write('\n')
elif (not rx_list):
    print("No reactions found!")
if wc_activate and (not wc_list):
    print("No Wannier center's movement found!")

#### HIGHLIGHTED TRAJECTORY FILE
if th_activate:                                         #-o h
    all_list=rx_list+wc_list
    if not all_list:
        print("Nothing to write!")
    else:
        xyz_highlight(filename, fname+"-H.xyz", all_list, "Si", "He")

#### BENCHMARKING
if benchmarking:                                        #-o b
    t3=timer()
    date1=datetime.now().strftime("%X %d/%m/%y")
    print("\n###################################################")
    print("##  VERSION          {0}".format(version))
    print("##  STARTED ON       {0}".format(date0))
    print("##  FINISHED ON      {0}".format(date1))
    print("##  TRAJECTORY FILE  {0}".format(filename))
    if args.radii:
        print("##  COVALENT RADII   {0}".format(cov_rad_file))
    else:
        print("##  COVALENT RADII   {0}".format('doi.org/10.1039/B801115J'))
    print("##  alpha            {0}".format(alpha))
    if wc_activate:
        print("##  beta             {0}".format(beta))
        print("##  gamma            {0}".format(gamma))
    print("##  REACTIONS        ELAPSED TIME  {0}".format(chrono(t1,t2)))
    if th_activate:
        print("##  H. TRAJECTORY    ELAPSED TIME  {0}".format(chrono(t2,t3)))
    print("##  GLOBAL           ELAPSED TIME  {0}".format(chrono(t0,t3)))
    print("###################################################")

