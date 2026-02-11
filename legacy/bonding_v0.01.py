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
dist_t0=timer()                                     #timer

traj_file=open(args.filename, 'r')
dist_file_name=fname+".DISTold.dat"
dist_file=open(dist_file_name, "w")
for s in range (traj_steps):
    dist_file.write(traj_file.readline())
    dist_file.write(traj_file.readline())
    in_step=list()
    for a in range (traj_N):
        in_step.append(traj_file.readline().split())
    out_dist=list()
    for i in range (traj_N):
        out_line=list()
        out_line.append(in_step[i][0])
        Xi=float(in_step[i][1])
        Yi=float(in_step[i][2])
        Zi=float(in_step[i][3])
        for j in range (traj_N):
            if (i>j):
                Dij=out_dist[j][i+1]
            elif (i==j):
                Dij=0
            elif (i<j):
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
traj_file.close()
dist_file.close()

dist_t1=timer()                                     #timer

#### CONNECTIVITY

cov_rad_data=open(cov_rad_file, "r")
cov_rad=list()
for l in cov_rad_data:
    line=l.split()
    if line[0]!="#":
        cov_rad.append(line)
cov_rad_data.close()

dist_file=open(dist_file_name, "r")
#bond_file_name=args.filename[:-4]+".BOND.dat"
#bond_file=open(bond_file_name, "w")
conn_file_name=fname+".CONNold.dat"
conn_file=open(conn_file_name, "w")

for s in range (traj_steps):
    conn_file.write(dist_file.readline())
    conn_file.write(dist_file.readline())
    in_dist=list()
    for a in range (traj_N):
        in_dist.append(dist_file.readline().split())
    out_conn=list()
    for i in range (traj_N):
        for element in cov_rad:
            if (in_dist[i][0]==element[1]):
                r_i=float(element[2])
        out_lcon="  "+in_dist[i][0].ljust(3," ")+str(i+1).rjust(6," ")
        for j in range (traj_N):
            for element in cov_rad:
                if (in_dist[j][0]==element[1]):
                    r_j=float(element[2])
            dist=float(in_dist[i][j+1])
            bdis=(r_i+r_j)*r_alpha
            if (dist<=bdis):
                bond=dist
            elif (dist>bdis):
                bond=-1
            if (bond>0):
                out_lcon+=str(j+1).rjust(6," ")
        conn_file.write(out_lcon+"\n")
dist_file.close()
#bond_file.close()
conn_file.close()

bond_t1=timer()                                     #timer

#### OUTPUT

conn_file=open(conn_file_name, "r")

report_name=fname+".REPORTold.txt"
report=open(report_name, "w")

#Step 0
conn_0=list()
conn_file.readline()
conn_file.readline()
for i in range (traj_N):
    conn_0.append(conn_file.readline().split())
    
#Steps 2 to the end

for s in range(1,traj_steps):
    conn_1=list()
    conn_file.readline()
    conn_file.readline()
#    cleav=list()
#    form=list()
    for i in range (traj_N):
        conn_1.append(conn_file.readline().split())
        if (conn_0[i]!=conn_1[i]):
            for bond in conn_0[i]:
                if (bond not in conn_1[i]) and (int(bond)>int(conn_0[i][1])):
                    AT1=conn_0[i][0]
                    ID1=conn_0[i][1]
                    AT2=conn_0[int(bond)-1][0]
                    ID2=conn_0[int(bond)-1][1]
                    report.write("Step "+str(s+1)+": Bond "+AT1+"("+ID1+")-"+AT2+"("+ID2+") is broken.\n")
            for bond in conn_1[i]:
                if (bond not in conn_0[i]) and (int(bond)>int(conn_0[i][1])):
                    AT1=conn_1[i][0]
                    ID1=conn_1[i][1]
                    AT2=conn_0[int(bond)-1][0]
                    ID2=conn_0[int(bond)-1][1]
                    report.write("Step "+str(s+1)+": Bond "+AT1+"("+ID1+")-"+AT2+"("+ID2+") is formed.\n")
    conn_0=conn_1
conn_file.close()


##########################################################
t1=timer()
date1=datetime.now().strftime("%X %d/%m/%y")
def chrono(start,end):
    return("{:.5f} minutes".format((end-start)/60).rjust(10," "))
report.write("\n############################################################\n")
report.write("##  STARTED ON       "+date0+"\n")
report.write("##  FINISHED ON      "+date1+"\n")
report.write("##  TRAJECTORY FILE  "+args.filename+"\n")
report.write("##  COVALENT RADII   "+cov_rad_file+"\n")
report.write("##  alpha            "+str(r_alpha)+"\n")
report.write("##  DISTANCES        ELAPSED TIME  "+chrono(dist_t0,dist_t1)+"\n")
report.write("##  CONNECTIVITIES   ELAPSED TIME  "+chrono(dist_t1,bond_t1)+"\n")
report.write("##  OUTPUT           ELAPSED TIME  "+chrono(bond_t1,t1)+"\n")
report.write("##  GLOBAL           ELAPSED TIME  "+chrono(t0,t1)+"\n")
report.write("############################################################\n")









