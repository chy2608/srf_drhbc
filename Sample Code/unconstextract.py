import os
import numpy as np
import matplotlib.pyplot as plt
import math
import shutil as shu

cdir = os.getcwd()
A = 126 #change

title = "                            neutron        proton         total\n"


target = [[-0.1, 0.15], [-0.15, 0.2]] #A-1, A+1

uncdef = []
E = []
initb = []

unconv = []
extracted = []
for t in target:
    for k in range(2):
        up = t[k]+0.15
        low = t[k]-0.1
        
        if up >= 0.65:
            up = 0.65
        else:
            pass

        if low <= -0.4:
            low = -0.4
        else:
            pass

        low = round(low, 2)
        up = round(up, 2)

        for i in np.arange(low, up, 0.05): #deformation values
            i = round(i,2)
            if math.isclose(i, 0) == True:
                i = 0
            Apath = os.path.join(cdir, "unconstrained", "newcodetest", f"{A}")
            if i in extracted:
                print("Data from " + os.path.join(Apath, f"b.{i}", "dir.out") + " already extracted to " 
                            + os.path.join(Apath, f"{A}u_bulk.txt"))
            else:
                with open(os.path.join(Apath, f"initb.{i}", "dir.out"), "r", encoding = 'utf8') as outfile:
                    data = outfile.readlines()
                    if data[-2] == " * No converged result in blocking calculations. *\n":
                        with open(os.path.join(Apath, f"{A}u_bulk.txt"), "a", encoding = 'utf8') as bulkfile:
                            bulkfile.write(f"initb = {i}\n")
                            bulkfile.write(" Did not converge\n")
                            bulkfile.write("-"*80+"\n")
                            unconv.append(i)
                            os.rename(os.path.join(Apath, f"initb.{i}"), os.path.join(Apath, f"initb.{i}u"))
                            bpath = os.path.join(Apath, f"initb.{i}")
                            os.mkdir(bpath)
                            shu.copy(os.path.join(Apath, f"initb.{i}u", "1drhbws"), bpath)
                            shu.copy(os.path.join(Apath, f"initb.{i}u", "dir_unconv.wel"), os.path.join(bpath, "dir.wel"))
                            src = os.path.join(Apath, f"initb.{i}u", "dir.dat")
                            with open(src, "r", encoding = 'utf8') as rf: #creates a new copy of dir.dat in each deformation subfolder with changed variables
                                with open(os.path.join(bpath, "dir.dat"), "w", encoding = 'utf8') as wf:
                                    for line in rf:
                                        if line == "inin     =    1                     ! initialization of potentials\n":
                                            wf.write("inin     =    0                     ! initialization of potentials\n")
                                        else:
                                            wf.write(line) 
                            os.chmod(os.path.join(bpath, "dir.dat"), 0o777)
                            print(f"initb={i} did not converge. Directory '{bpath}' created. All filed succcesfully appended with permissions -rwxrwxrwx")
                    else:
                        initb.append(i)
                        ind = len(data) - 1 - data[::-1].index(title)
                        with open(os.path.join(Apath, f"{A}u_bulk.txt"), "a", encoding = 'utf8') as bulkfile:
                            bulkfile.write(f"initb = {i}\n")
                            bulkfile.write(title + "\n")
                            bulkfile.write(data[ind+9] + "\n")
                            bulkfile.write(data[ind+22]+ "\n")
                            E.append(float(data[ind+22][55:].rstrip()))
                            uncdef.append(float(data[ind+9][58:].rstrip()))
                            bulkfile.write("-"*80+"\n")
                            print("Data from " + os.path.join(Apath, f"b.{i}", "dir.out") + " appended to " 
                                + os.path.join(Apath, f"{A}u_bulk.txt"))
                extracted.append(i)

print("")
minE = min(E)
minb = uncdef[E.index(minE)]

if len(unconv) != 0:
    for b in unconv:
        with open(os.path.join(cdir, "unconstrained", "shells2", f"EU{A}unconv.sh"), "a") as shell:
                shell.write(f"(cd {os.path.join(Apath, f'initb.{b}')} && ./1drhbws)\n") #creates shell executable which runs all unconstrained deformation 1drhbws 
    os.chmod(os.path.join(cdir, "unconstrained", "shells2", f"EU{A}unconv.sh"), 0o777) 
    print(f"Shell file for unconverged initb created at '{os.path.join(cdir, 'unconstrained', 'shells2', f'EU{A}unconv.sh')}' with permissions -rwxrwxrwx")  

with open(os.path.join(Apath, f"{A}u_bulk.txt"), "a", encoding = 'utf8') as bulkfile:
    bulkfile.write("Minimum E, b: \n{:.6f},{}".format(minE, minb))

gs = initb[E.index(min(E))]


with open(os.path.join(Apath, f"initb.{gs}", "dir.out"), "r", encoding = 'utf8') as gsoutfile:
    data = gsoutfile.readlines()
    ind = len(data) - 1 - data[::-1].index(title)
    with open(os.path.join(cdir, "groundstates2.txt"), "a", encoding = 'utf8') as gsdata:
        gsdata.write("{:<8.0f}{:<8.0f}{:<10.2f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.2f}{:<10.2f}{:<10}{}\n".format(
            A #A
            , int(A-63) #N 
            , abs(float(data[ind+22][55:].rstrip())) #Eb
            , abs(float(data[ind+23][58:].rstrip())) #Erot
            , float(data[ind+6][29:37]) #Rn
            , float(data[ind+6][44:52]) #Rp
            , float(data[ind+6][59:].rstrip()) #Rm
            , float(data[ind+7][44:52]) #Rc
            , float(data[ind+9][28:37].strip()) #bn2
            , float(data[ind+9][43:52].strip()) #bp2
            , float(data[ind+9][58:].strip()) #b2
            , float(data[ind+2][25:37].strip()) #ln
            , float(data[ind+2][40:].strip()) #lp
            , str(data[ind+26][-10:].strip()) #mpi
            , "          "))
        print("Data from " + os.path.join(Apath, f"initb.{gs}", "dir.out") + " appended to " 
                + os.path.join(cdir, "groundstates2.txt"))
print("")

with open(os.path.join(cdir, "bulk2.txt"), "a", encoding = 'utf8') as bulkdata:
        bulkdata.write("{:<8}{:<8}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format(
            "A", "N", "init b2", "Etot+cm", "Erot", "Rn", "Rp", "Rm", "Rc", "bn2", "bp2", "b2", "ln", "lp", "mpi"))

bc = 0
for b in sorted([float(k) for k in initb]):
    if b == 0:
        b = int(b)
    with open(os.path.join(Apath, f"initb.{b}", "dir.out"), "r", encoding = 'utf8') as boutfile:
        if b == gs:
            bstr = str(b)+"*"
        else:
            bstr = str(b)
        data = boutfile.readlines()
        ind = len(data) - 1 - data[::-1].index(title)
        with open(os.path.join(cdir, "bulk2.txt"), "a", encoding = 'utf8') as bulkdata:
            if bc == 0:
                bulkdata.write("{:<8.0f}{:<8.0f}{:<10}{:<10.2f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.2f}{:<10.2f}{:<10}\n".format(
                    A #A
                    , int(A-63) #N 
                    , str(bstr) # initb
                    , float(data[ind+22][55:].rstrip()) #Etotcm
                    , float(data[ind+23][58:].rstrip()) #Erot
                    , float(data[ind+6][29:37]) #Rn
                    , float(data[ind+6][44:52]) #Rp
                    , float(data[ind+6][59:].rstrip()) #Rm
                    , float(data[ind+7][44:52]) #Rc
                    , float(data[ind+9][28:37].strip()) #bn2
                    , float(data[ind+9][43:52].strip()) #bp2
                    , float(data[ind+9][58:].strip()) #b2
                    , float(data[ind+2][25:37].strip()) #ln
                    , float(data[ind+2][40:].strip()) #lp
                    , str(data[ind+26][-10:].strip()) #mpi
                    ))
            else:
                bulkdata.write("{:<8}{:<8}{:<10}{:<10.2f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.3f}{:<10.2f}{:<10.2f}{:<10}\n".format(
                    " "
                    , " " 
                    , str(bstr) # initb
                    , float(data[ind+22][55:].rstrip()) #Etotcm
                    , float(data[ind+23][58:].rstrip()) #Erot
                    , float(data[ind+6][29:37]) #Rn
                    , float(data[ind+6][44:52]) #Rp
                    , float(data[ind+6][59:].rstrip()) #Rm
                    , float(data[ind+7][44:52]) #Rc
                    , float(data[ind+9][28:37].strip()) #bn2
                    , float(data[ind+9][43:52].strip()) #bp2
                    , float(data[ind+9][58:].strip()) #b2
                    , float(data[ind+2][25:37].strip()) #ln
                    , float(data[ind+2][40:].strip()) #lp
                    , str(data[ind+26][-10:].strip()) #mpi
                    ))
            
            print("Data from " + os.path.join(Apath, f"initb.{b}", "dir.out") + " appended to " 
                    + os.path.join(cdir, "bulk2.txt"))
            bc += 1    
with open(os.path.join(cdir, "bulk2.txt"), "a", encoding = 'utf8') as bulkdata:
    bulkdata.write("-"*141+"\n")

print("")

plt.plot(uncdef, E, "v", label = "Unconstrained")

plt.legend(loc = "best")
plt.xlabel(r"Deformation ($\beta$)")
plt.ylabel(r"Etot + Ecm ($\mathrm{MeV}$)")
plt.title("Total energy against quadrupole deformation ((Un)constrained)")
plt.savefig(os.path.join(Apath, f"{A}_unconstrained.jpg"), dpi = 2000)
print(f"Plot generated at '{os.path.join(Apath, f'{A}_unconstrained.jpg')}'")
print("")
print("End")