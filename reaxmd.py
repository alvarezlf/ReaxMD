#### ReaxMD
#### Covalent radii by B. Cordero et. al., Dalton Transactions, 2008, 21. https://doi.org/10.1039/B801115J

#### FUNCTIONS
def chrono(start,end):
        o_timer=end-start
        if o_timer>=60:
            return("{:.3f} minutes".format(o_timer/60))
        else:
            return("{:.3f} seconds".format(o_timer))

def step2coor(file, N):                     #file, 1 step -> coord. matrix
    file.readline()
    file.readline()
    coor=list()
    for i in range(1,N+1):
        line=file.readline().split()
        line.insert(0,i)
        coor.append(line)
    return coor
#[[ID1, "A1", "X1", "Y1", "Z1"],...]

def strider(file, s_val):
    if s_val > 0 and type(s_val)==int:
        for i in range ((s_val-1)*(traj_N+2)):
            file.readline()
    else:
        print('Error: wrong stride value.')
        sys.exit()
    return

#### BONDING FUNCTION
#coord matrix -> atoms close to Wanniers
def wannier(coor_w, coor_a, cutoff):    #vectorized
    coords_w = coor_w[:, 2:5].astype(float)
    coords_a = coor_a[:, 2:5].astype(float)
    delta = (coords_a[np.newaxis, :, :] - coords_w[:, np.newaxis, :] + box_half) % box_size - box_half
    dist_matrix = np.sqrt(np.sum(delta**2, axis=2))
    conn_matrix = dist_matrix <= cutoff[np.newaxis, :]
    wann_step = [[tuple(coor_w[i, :2])] + [tuple(coor_a[j, :2]) for j in np.where(conn_matrix[i])[0]] for i in range(len(coor_w))]
    return wann_step
#[[ID1, "A1", d11, d12, d13...],[ID2, "A2",d22, d23...]...]

def coor2bond(coor_step, N, cutoff):    #vectorized
    coords = np.array([i[2:5] for i in coor_step]).astype(float)
    delta = (coords[:, None, :] - coords[None, :, :] + box_half) % box_size - box_half
    dist_matrix = np.sqrt(np.sum(delta**2, axis=2))
    bond_matrix = (dist_matrix <= cutoff) & (np.triu(np.ones((N, N)), k=1).astype(bool))
    atom_list=[tuple(i[:2]) for i in coor_step]
    bond_step = [[atom_list[i]] + [tuple(atom_list[j]) for j in np.where(bond_matrix[i])[0]] for i in range(N)]
    return bond_step

def bond2elem(c_elem, c_rest, cutoff):
    coor_elem = c_elem[:, 2:5].astype(float)
    coor_rest = c_rest[:, 2:5].astype(float)
    delta = (coor_rest[np.newaxis, :, :] - coor_elem[:, np.newaxis, :] + box_half) % box_size - box_half
    dist_matrix = np.sqrt(np.sum(delta**2, axis=2))
    bond_matrix = dist_matrix <= cutoff
    bond_step = [[tuple(c_elem[i, :2])] + [tuple(c_rest[j, :2]) for j in np.where(bond_matrix[i])[0]] for i in range(len(c_elem))]
    return bond_step

#[[ID1, A1, [a1, id1], [a2, id2]...]...]

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
                    reac.append(list([ID1, AT1, ID2, AT2, -1]))
            for bond in b2[i]:
                if (bond not in b1[i]):
                    AT1=b2[i][0][1]
                    ID1=b2[i][0][0]
                    AT2=bond[1]
                    ID2=bond[0]
                    reac.append(list([ID1, AT1, ID2, AT2, +1]))
    return(reac)

#[[Step, n1, a1, n2, a2, "formed/broken"],...]

def reax_filter(reac, threshold = 100):
    discart = set()
    output = []
    N = len(reac)
    for i in range(N):
        for j in range(i+1, N):
            if (tuple(reac[i]) not in discart) and (reac[j][0] - reac[i][0] <= threshold) and (reac[i][1:5] == reac[j][1:5]) and (reac[i][5] + reac[j][5] == 0):
                discart.add (tuple (reac[i]))
                discart.add (tuple(reac[j]))
                break
        if (tuple(reac[i]) not in discart):
            output.append(reac[i])
    return output

def print_reac(reac):
    print("step {0}{6}: bond ({1}){2}-{4}({3}) {5}".format(reac[0], reac[1], reac[2], reac[3], reac[4], 'formed' if reac[5]==1 else 'broken', ' at {0:.3f} ps'.format(reac[0]*timegap) if args.timestep else ''))
    return

def print_wann(wann):
    print("step {0}{6}: Wannier center {1} {5} {4}({3})".format(wann[0], wann[1], wann[2], wann[3], wann[4], 'approached' if wann[5]==1 else 'separated from', ' at {0:.3f} ps'.format(wann[0]*timegap) if args.timestep else ''))
    return

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
        print("\nError: Dummy atom \"X\" found in {0}!".format(filename))
        sys.exit()
    if wc_activate and ("X" not in atoms):
        print("\nError: No Wannier centers found in {0}!".format(filename))
        sys.exit()

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

def reactive_molecules(list_r, list_rx):
    reac_atoms=set([tuple([i[1], i[2]]) for i in list_rx] + [tuple([i[3], i[4]]) for i in list_rx])
    set_rx=set()
    for i in reac_atoms:
        for j in list_r:
            if i in j:
                set_rx.add(tuple(j))
                break
    list_rx_unsorted=list(set_rx)
    list_rx=[]
    for i in range(1,traj_N+1):
        for m in list_rx_unsorted:
            if m[0][0]==i:
                list_rx.append(m)
                break
    return(list_rx)

def net_reaction(mole_0, mole_1, rx_list):
    if not rx_list:
        return('\nNo reaction found!')
    else:
        form_0=mol_formula(mole_0)
        form_1=mol_formula(mole_1)
        for i in form_0:
            if i not in form_1:
                form_1[i]=0
        for i in form_1:
            if i not in form_0:
                form_0[i]=0
        form_Rx={}
        Reac=''
        Prod=''
        for i in form_1:
            form_Rx[i]=form_1[i]-form_0[i]
        for i in form_Rx:
            if form_Rx[i]>0:
                Prod+=' + {0} {1}'.format(form_Rx[i], i)
                Prod=Prod.replace('1 ','')
            elif form_Rx[i]<0:
                Reac+='{0} {1} + '.format(-form_Rx[i], i)
                Reac=Reac.replace('1 ','')
        net_rx=[i for i in form_Rx if form_Rx[i]!=0]
        if net_rx:
            return('\nNet reaction:\n{0} -> {1}'.format(Reac[:-2], Prod[2:]))
        else:
            return('\nNo net reaction found!')

def reaxmd():
    rx_list=list()
    wc_list=list()
#### ANALYSIS
    print('trajectory file: {0}'.format(filename))
    print('cell sides: {0} {1} {2} {3}'.format(args.A, args.B, args.C, 'a.u.' if args.atomic_units else 'Å'))
    print('atoms: {0}'.format(traj_N))
    print('steps: {0}'.format(traj_steps))
    if args.stride:
        print ("stride: {0}".format(args.stride))
    if args.timestep:
        print("time step: {0} a.u.".format(args.timestep))

#### MAIN
    if args.element and not args.reaction_only:                                        #-e
        with open(filename, 'r') as traj_file:
            s=1
            coor_0=step2coor(traj_file,traj_N)
            checkdummy(coor_0)
            if len(ELEM)>2 or not ELEM:
                print("\nPlease specify 1 or 2 elements as \"-e E0\" or \"-e E1-E2\"")
                sys.exit()
            elif len(ELEM)==2:
                coor_0a1=np.array([i for i in coor_0 if i[1] == ELEM[0]])
                coor_0a2=np.array([i for i in coor_0 if i[1] == ELEM[1]])
                types_elem = coor_0a1[:, 1]
                types_rest = coor_0a2[:, 1]
                r_elem = np.array([cov_rad[t] for t in types_elem])
                r_rest = np.array([cov_rad[t] for t in types_rest])
                cutoff_matrix = (r_elem[:, np.newaxis] + r_rest) * alpha
                bond_0=bond2elem(coor_0a1, coor_0a2, cutoff_matrix)
            elif len(ELEM)==1:
                coor_0e=np.array([i for i in coor_0 if i[1] in ELEM])
                coor_0r=np.array([i for i in coor_0 if i[1] not in ELEM])
                types_elem = coor_0e[:, 1]
                types_rest = coor_0r[:, 1]
                r_elem = np.array([cov_rad[t] for t in types_elem])
                r_rest = np.array([cov_rad[t] for t in types_rest])
                cutoff_matrix = (r_elem[:, np.newaxis] + r_rest) * alpha
                bond_0=bond2elem(coor_0e, coor_0r, cutoff_matrix)
            if args.stride:
                strider(traj_file, args.stride)
                s+=args.stride
            else:
                s+=1
            while s<=traj_steps:
#            for s in range (2,traj_steps+1):
                coor_1=step2coor(traj_file,traj_N)
                if ELEM and len(ELEM)==2:
                    coor_1a1=np.array([i for i in coor_1 if i[1]==ELEM[0]])
                    coor_1a2=np.array([i for i in coor_1 if i[1]==ELEM[1]])
                    bond_1=bond2elem(coor_1a1, coor_1a2, cutoff_matrix)
                    rx_1=bond2reac(bond_0, bond_1)
                    rx=[]
                    for i in rx_1:
                        if [i[2],i[3],i[0],i[1],i[4]] not in rx:
                            rx.append(i)
                elif ELEM and len(ELEM)==1:
                    coor_1e=np.array([i for i in coor_1 if i[1] in ELEM])
                    coor_1r=np.array([i for i in coor_1 if i[1] not in ELEM])
                    bond_1=bond2elem(coor_1e, coor_1r, cutoff_matrix)
                    rx=bond2reac(bond_0, bond_1)
                if rx:
                    if not rx_list:
                        print('')
                    for i in rx:
                        i.insert(0,s)
                        if not out_filter:
                            print_reac(i)
                        rx_list.append(i)
                coor_0=coor_1
                bond_0=bond_1
                sys.stdout.flush()
                if args.stride:
                    strider(traj_file, args.stride)
                    s+=args.stride
                else:
                    s+=1

    if wc_activate and not args.reaction_only:                                       #-o w
        with open(filename, 'r') as traj_file:
            s=1
            coor_0=step2coor(traj_file,traj_N)
            checkdummy(coor_0)
            coor_w0=np.array([i for i in coor_0 if i[1]=="X"])
            coor_0=np.array([i for i in coor_0 if i[1]!="X"])
            atom_types = coor_0[:, 1]
            covrad_matrix = np.array([cov_rad[t] for t in atom_types])
            cutoff_matrix = (covrad_matrix[:, None] + covrad_matrix) * alpha
            cutoff_wannier = beta * covrad_matrix + gamma
            N_at=len(coor_0)
            bond_0=coor2bond(coor_0, N_at, cutoff_matrix)
            wann_0=wannier(coor_w0, coor_0, cutoff_wannier)
            if args.stride:
                strider(traj_file, args.stride)
                s+=args.stride
            else:
                s+=1
            while s<=traj_steps:
#            for s in range (2,traj_steps+1):
                coor_1=step2coor(traj_file,traj_N)
                coor_w1=np.array([i for i in coor_1 if i[1]=="X"])
                coor_1=np.array([i for i in coor_1 if i[1]!="X"])
                bond_1=coor2bond(coor_1, N_at, cutoff_matrix)
                wann_1=wannier(coor_w1, coor_1, cutoff_wannier)
                rx=bond2reac(bond_0, bond_1)
                if rx:
                    if not rx_list:
                        print('')
                    for i in rx:
                        i.insert(0,s)
                        if not out_filter:
                            print_reac(i)
                        rx_list.append(i)
                wann=bond2reac(wann_0, wann_1)

                if wann:
                    for i in wann:
                        i.insert(0,s)
                        wc_list.append(i)
                coor_0=coor_1
                bond_0=bond_1
                coor_w0=coor_w1
                wann_0=wann_1
                sys.stdout.flush()
                if args.stride:
                    strider(traj_file, args.stride)
                    s+=args.stride
                else:
                    s+=1

        if not wc_list:
            print("\nNo Wannier center's movement found!")
        elif not out_filter:
            print()
            for i in wc_list:
                print_wann(i)
                    
    if normal:
        with open(filename, 'r') as traj_file:
            s=1
            coor_0=step2coor(traj_file, traj_N)
            checkdummy(coor_0)
            atom_types = [i[1] for i in coor_0]
            covrad_matrix = np.array([cov_rad[t] for t in atom_types])
            cutoff_matrix = (covrad_matrix[:, None] + covrad_matrix) * alpha
            bond_0=coor2bond(coor_0, traj_N, cutoff_matrix)
            if args.stride:
                strider(traj_file, args.stride)
                s+=args.stride
            else:
                s+=1
            while s<=traj_steps:
            #for s in range (2,traj_steps+1):
                coor_1=step2coor(traj_file, traj_N)
                bond_1=coor2bond(coor_1, traj_N, cutoff_matrix)
                rx=bond2reac(bond_0,bond_1)
                if rx:
                    if not rx_list:
                        print('')
                    for i in rx:
                        i.insert(0,s)
                        if not out_filter:
                            print_reac(i)
                        rx_list.append(i)
                coor_0=coor_1
                bond_0=bond_1
                sys.stdout.flush()
                if args.stride:
                    strider(traj_file, args.stride)
                    s+=args.stride
                else:
                    s+=1


    if out_filter:
        if rx_list:
            rx_list = reax_filter (rx_list)
            for i in rx_list:
                print_reac(i)
        if wc_list:
            rx_atom=set([i[1] for i in rx_list]+[i[3] for i in rx_list])
            wc_f=[w for w in wc_list if w[3] in rx_atom]
            wc_list = reax_filter(wc_f)
            print()
            for i in wc_list:
                    print_wann(i)

#### NET REACTION
    with open(filename, 'r') as traj_file:
        coor_0=step2coor(traj_file, traj_N)
        checkdummy(coor_0)
        wann_found = 'X' in [i[1] for i in coor_0]
        if wann_found:
            coor_nw=[i for i in coor_0 if i[1]!='X']
            coor_0=coor_nw
        atom_types = [i[1] for i in  coor_0]
        covrad_matrix = np.array([cov_rad[t] for t in atom_types])
        cutoff_matrix = (covrad_matrix[:, None] + covrad_matrix) * alpha
        bond_0=coor2bond(coor_0, len(coor_0), cutoff_matrix)
        strider(traj_file, traj_steps-1)
        coor_1=step2coor(traj_file, traj_N)
        if wann_found:
            coor_nw=[i for i in coor_1 if i[1]!='X']
            coor_1=coor_nw
        bond_1=coor2bond(coor_1, len(coor_1), cutoff_matrix)
        rx=bond2reac(bond_0,bond_1)
        mole_0=join_molecules(bond_0)
        mole_1=join_molecules(bond_1)
        ####
        if rx and args.reaction_only:
            for i in rx:
                i.insert(0,traj_steps)
                rx_list.append(i)
        ####
    print('{0}'.format(net_reaction(mole_0, mole_1, rx)))
    
#### WRITE MOLECULES FILE
    if args.options and ('m' in write_opt) and rx_list:     #-o m
        react_list=reactive_molecules(mole_0, rx_list)
        with open(fname+'-molecules.txt', 'w') as mole_file:
            mole_file.write('{0}\n{1}\n'.format(filename, traj_N))
            mole_file.write(r'# List of reactive molecules generated automatically by ReaxMD.'+'\n')
            if out_filter or args.reaction_only:
                mole_file.write(r'# Warning: List created by reduced analysis and may not represent all reactive molecules.'+'\n')
            for i in react_list:
                for j in i:
                    mole_file.write(str(j[0])+'\n')
                mole_file.write('\n')

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
        print("##  ReaxMD Version   {0}".format(version_val))
        if args.radii:
            print("##  COVALENT RADII   {0}".format(cov_rad_file))
        else:
            print("##  COVALENT RADII   {0}".format('doi.org/10.1039/B801115J'))
        print("##  ALPHA VALUE      {0}".format(alpha))
        if wc_activate:
            print("##  BETA VALUE       {0}".format(beta))
            print("##  GAMMA VALUE      {0}".format(gamma))
        if out_filter:
            print("##  FAST OSCILATIONS FILTERED")
        print("##  EXECUTION TIME   {0}".format(chrono(t0,t3)))
        print("###################################################")

    return rx_list, wc_list

if __name__ == "__main__":
    version_val="260107"

    import os
    import sys
    import argparse
    import numpy as np
    from collections import Counter
#### INPUT
    parser=argparse.ArgumentParser()
    parser.add_argument("trajectory_files", nargs = '*', help="Path(s) to trajectory file(s).")
    parser.add_argument("A", type=float, help="A length of the cell in Angstrom.")
    parser.add_argument("B", type=float, help="B length of the cell in Angstrom.")
    parser.add_argument("C", type=float, help="C length of the cell in Angstrom.")
    parser.add_argument("-au", "--atomic_units", action="store_true", help="Specifies that the cell side lengths are in atomic units.")
    parser.add_argument("-a", "--alpha", type=float, default=1.30, help="Alpha value for bond lengths as a*(r_i + r_j). Default value: 1.30")
    parser.add_argument("-b", "--beta", type=float, default=1.30, help="Beta value for Wannier centers proximities as (b*r_i + g). Default value: 1.30")
    parser.add_argument("-g", "--gamma", type=float, default=0.40, help="Gamma value for Wannier centers proximities as (b*r_i + g). Default value: 0.40")
    parser.add_argument("-e", "--element", type=str, help=("\"-e El\" analyzes only the bonds involving element El.\n"
                        "\"-e E1-E2\" analyzes only bonds between elements E1 and E2."))
    parser.add_argument("-o", "--options", type=str, help=("\"-o hmwfi\", where: h writes a trajectory file \"trajectory_file-H.xyz\", where the atoms involved in reactions are replaced by Si and the Wannier centers (include \"w\" as option) by He. m writes a file \"trajectory-molecules.txt\" with the reactive molecules. w analyzes the Wannier centers. f filters equilibrium reactions within 100 MD steps. i prints the working parameters information at the end."))
    parser.add_argument("-r", "--radii", type=str, help="Alternative full path to covalent radii data.")
    parser.add_argument("-ro", "--reaction_only", action='store_true', help="Print reaction only.")
    parser.add_argument("-s", "--stride", type=int, help="Jumps the specified amount of steps.")
    parser.add_argument("-t", "--timestep", type=float, help="Prints the time in ps next to the step number. Timestep = (MD timestep) x (MD steps between trajectory steps)")
    parser.add_argument('-v', '--version', action='version', version='ReaxMD version {0}'.format(version_val), help='Print version information')
    args=parser.parse_args()

    #### CONFIGURATION
    box_size = np.array([args.A, args.B, args.C])
    if args.atomic_units:
        box_size = box_size * 0.529177210903
    box_half = box_size/2

    alpha=args.alpha
    beta=args.beta
    gamma=args.gamma

    write_opt=args.options

    benchmarking=(args.options) and ("i" in write_opt)      #timer
    if benchmarking:
        from datetime import datetime
        from timeit import default_timer as timer
        t0=timer()
        date0=datetime.now().strftime("%X %d/%m/%y")

    wc_activate=(args.options) and ("w" in write_opt)
    out_filter = write_opt and ("f" in write_opt)
    th_activate=(args.options) and ("h" in write_opt)
    normal= (not args.element) and (not wc_activate) and (not args.reaction_only)

    if args.element:
        ELEM=args.element.split("-")
        if wc_activate:
            print('Error: Wannier analysis is incompatible with element filtering!')
            sys.exit()

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

    if args.timestep:
        timegap=args.timestep*2.41888e-5

    for filename in args.trajectory_files:
        fname = filename.replace(".xyz","")
        with open(filename, 'r') as traj_file:
            traj_N = int(traj_file.readline().strip())
            traj_lines = sum(1 for i in traj_file)+1
        traj_steps = int(traj_lines/(traj_N + 2))

        if traj_steps==1:
            print('Error: Only 1 step found in {0}'.format(filename))
            sys.exit()
        try:
            reactions, wanniers = reaxmd()
        except KeyboardInterrupt:
            print("Keyboard Interrupt!")
    
#### Have a nice day
####  ~ L. Álvarez ~
