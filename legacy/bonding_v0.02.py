#### Timer
from datetime import datetime
from timeit import default_timer as timer
t0=timer()                                          #timer
date0=datetime.now().strftime("%X %d/%m/%y")

#### Enter parameters
cov_rad_file="/home/alvarez/bin/covalent_radii.ssv" ### Absolute path to the covalent radii database
r_alpha=1.3                                         ### Margin for comparing the bond distances vs r_cov

#### Collect data
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
traj_lines=1
traj_N=int(traj_file.readline().replace("\n", ""))
for i in traj_file:
    traj_lines+=1
traj_steps=traj_lines/int(traj_N+2)
traj_file.close()

#### DISTANCES
t_dist_0=timer()                                     #timer

traj_file=open(args.filename, 'r')
dist_file_name=fname+".DISTv2.dat"
dist_file=open(dist_file_name, "w")
in_step=list()
for s in range (traj_steps):
    dist_file.write(traj_file.readline())
    dist_file.write(traj_file.readline())
    for a in range (traj_N):
        in_step.append(traj_file.readline().split())
    out_dist=list()
    for i in range (traj_N):
        out_line=list()
        out_line.append(in_step[i][0])
        Xi=float(in_step[i][1])
        Yi=float(in_step[i][2])
        Zi=float(in_step[i][3])
        for j in range (i,traj_N):
            DXj=float(in_step[j][1])-Xi
            DYj=float(in_step[j][2])-Yi
            DZj=float(in_step[j][3])-Zi
            while DXj>A/2:
                DXj-=A
            while DXj<-A/2:
                DXj+=A
            while DYj>B/2:
                DYj-=B
            while DYj<-B/2:
                DYj+=B
            while DZj>C/2:
                DZj-=C
            while DZj<-C/2:
                DZj+=C
            Dij=sqrt(DXj**2+DYj**2+DZj**2)
            out_line.append(round(Dij,3))
        out_dist.append(out_line)
    for l in out_dist:
        txt='  '+l[0].ljust(3," ")
        for i in l[1:]:
            txt+="{:.3f}".format(i).rjust(10," ")
        dist_file.write(txt+"\n")
    in_step=list()
traj_file.close()
dist_file.close()

t_dist_1=timer()                                     #timer

#### CONNECTIVITY

cov_rad_data=open(cov_rad_file, "r")
cov_rad=list()
for l in cov_rad_data:
    line=l.split()
    if line[0]!="#":
        cov_rad.append(line)
cov_rad_data.close()

dist_file=open(dist_file_name, "r")
conn_file_name=fname+".CONNv2.dat"
conn_file=open(conn_file_name, "w")
in_dist=list()
for s in range (traj_steps):
    conn_file.write(dist_file.readline())
    conn_file.write(dist_file.readline())
    for a in range (traj_N):
        in_dist.append(dist_file.readline().split())
    out_conn=list()
    for i in range (traj_N):
        for element in cov_rad:
            if (in_dist[i][0]==element[1]):
                r_i=float(element[2])
        out_lcon="  "+in_dist[i][0].ljust(3," ")+str(i+1).rjust(6," ")
        for j in range (traj_N-i):
            for element in cov_rad:
                if (in_dist[j+i][0]==element[1]):
                    r_j=float(element[2])
            dist=float(in_dist[i][j+1])
            bdis=(r_i+r_j)*r_alpha
            if (dist<=bdis):
                bond=dist
            elif (dist>bdis):
                bond=-1
            if (bond>0):
                out_lcon+=str(i+j+1).rjust(6," ")
        conn_file.write(out_lcon+"\n")
    in_dist=list()
dist_file.close()
conn_file.close()

t_bond=timer()                                     #timer

#### OUTPUT

conn_file=open(conn_file_name, "r")

report_name=fname+".REPORTv2.txt"
report=open(report_name, "w")

#Step 0
conn_0=list()
conn_file.readline()
conn_file.readline()
for i in range (traj_N):
    conn_0.append(conn_file.readline().split())
    
#Steps 2 to the end
rx_list=list()

for s in range(1,traj_steps):
    conn_1=list()
    conn_file.readline()
    conn_file.readline()
    for i in range (traj_N):
        conn_1.append(conn_file.readline().split())
        if (conn_0[i]!=conn_1[i]):
            for bond in conn_0[i]:
                if (bond not in conn_1[i]):
                    AT1=conn_0[i][0]
                    ID1=int(conn_0[i][1])
                    AT2=conn_0[int(bond)-1][0]
                    ID2=int(conn_0[int(bond)-1][1])
                    rx_list.append(tuple([s+1, ID1, AT1, ID2, AT2, "broken"]))
                    report.write("Step {0}: Bond {2}({1})-{4}({3}) is {5}\n".format(rx_list[-1][0], rx_list[-1][1], rx_list[-1][2], rx_list[-1][3], rx_list[-1][4], rx_list[-1][5]))
            for bond in conn_1[i]:
                if (bond not in conn_0[i]):
                    AT1=conn_1[i][0]
                    ID1=int(conn_1[i][1])
                    AT2=conn_0[int(bond)-1][0]
                    ID2=int(conn_0[int(bond)-1][1])
                    rx_list.append(tuple([s+1, ID1, AT1, ID2, AT2, "formed"]))
                    report.write("Step {0}: Bond {2}({1})-{4}({3}) is {5}.\n".format(rx_list[-1][0], rx_list[-1][1], rx_list[-1][2], rx_list[-1][3], rx_list[-1][4], rx_list[-1][5]))
    conn_0=conn_1
conn_file.close()
report.close()
if rx_list==[]:
    report.write("No reactions found!\n")
t_outp=timer()                                     #timer

#### HIGHLIGHTED TRAJECTORY FILE
if rx_list!=[]:
    atom_list=list()
    for r in rx_list:
        atom_list.append([r[0],r[1],r[2]])
        atom_list.append([r[0],r[3],r[4]])
    atom_list.sort()
    traj_file=open(args.filename, "r")
    htra_file_name=fname+"-HTRAv2.xyz"
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

t_htra=timer()
##########################################################
t1=timer()
date1=datetime.now().strftime("%X %d/%m/%y")
def chrono(start,end):
    o_timer=end-start
    if o_timer>=60:
        return("{:.3f} minutes".format(o_timer/60))
    else:
        return("{:.3f} seconds".format(o_timer))

report=open(report_name, "a")
report.write("\n############################################################\n")
report.write("##  STARTED ON       {0}\n".format(date0))
report.write("##  FINISHED ON      {0}\n".format(date1))
report.write("##  TRAJECTORY FILE  {0}\n".format(args.filename))
report.write("##  COVALENT RADII   {0}\n".format(cov_rad_file))
report.write("##  alpha            {0}\n".format(r_alpha))
report.write("##  DISTANCES        ELAPSED TIME  {0}\n".format(chrono(t_dist_0,t_dist_1)))
report.write("##  CONNECTIVITIES   ELAPSED TIME  {0}\n".format(chrono(t_dist_1,t_bond)))
report.write("##  OUTPUT           ELAPSED TIME  {0}\n".format(chrono(t_bond,t_outp)))
report.write("##  H. TRAJECTORY    ELAPSED TIME  {0}\n".format(chrono(t_outp,t_htra)))
report.write("##  GLOBAL           ELAPSED TIME  {0}\n".format(chrono(t0,t1)))
report.write("############################################################\n")
report.close()






