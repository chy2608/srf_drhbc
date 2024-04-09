import numpy as np
import os
import matplotlib.pyplot as plt
import time

st = time.time()

ferm = 2 # 1=neutron 2=proton
component = False
dim2 = False
param = True # if True, dim2 must be True
trace = 2 # 0 = no trace, 1 = trace, 2 = both
Alist = [197] # if component == True or dim2 == True, only components of A=Alist[0] will be plotted

x = []
pden = []
nden = []

for A in Alist:
    Afile = os.path.join(os.getcwd(), "unconstrained", "data", "EU", f"{A}")
    #Afile = os.path.join(os.getcwd(), "unconstrained", "data", f"{A}")
    #searchFile = os.path.join(Afile, f"{A}u_bulk.txt")
    #with open(searchFile, "r") as searchData:
    #    data = searchData.readlines()
    #    target = f" Etot + Ecm                                            {data[-1][:12]}\n"
    #    for line in data:
    #        if line == target:
    #            minB = round(float(data[data.index(line)-5][7:].strip()), 2)
    #            if minB == 0:
    #                minB = 0

    title = "      x (fm)       ROS n       ROV n       ROS p       ROV p\n"

    xA = []
    ndenA = []
    pdenA = []

    #with open(os.path.join(Afile, f"initb.{minB}", "dir.out")) as outfile:
    with open(os.path.join(Afile, f"gs", "dir.out")) as outfile:
        data = outfile.readlines()
        start = data.index(title)+1
        end = start+201
        l = 0
        while l < 4:
            dens = data[start:end]
            xl = []
            pdenl = []
            ndenl = []
            for line in dens:
                if line[2:12] == "0.1000E-03":
                    l += 1
                xl.append(float(line[:8].strip())*10**(float(line[9:13].strip())))
                ndenl.append(float(line[13:20].strip())*10**(float(line[21:25].strip())))
                pdenl.append(float(line[37:44].strip())*10**(float(line[45:49].strip())))
            xA.append(xl)
            ndenA.append(ndenl)
            pdenA.append(pdenl)
            endt = end
            start = end
            end = endt+201
    
    x.append(xA)
    nden.append(ndenA)
    pden.append(pdenA)
                

if dim2 == False:
    fig = plt.figure(figsize=(7, 10))

    if component == True:
        for i in range(4):
            if ferm == 1:
                new = nden[0][i]
                plt.plot(x[0][i], new, label = r"$\lambda$ = " + f"{i*2}", lw = 0.7)
            else:
                new = pden[0][i]
                plt.plot(x[0][i], new, label = r"$\lambda$ = " + f"{i*2}", lw = 0.7)
    else:
        xrange = [i for i in x[0][0] if i < 0.5]
        
        if ferm == 1:
            sumN = []
            for A in nden:
                cur = []
                for i in range(201):
                    suml = []
                    for l in range(4):
                        suml.append(A[l][i])
                    cur.append(sum(suml))
                sumN.append(cur)
            for k in Alist:
                plt.plot(x[0][0], sumN[list(Alist).index(k)], "--", label = f"{k}", lw = 0.7)
        else:
            sumP = []
            for A in pden:
                cur = []
                for i in range(201):
                    suml = []
                    for l in range(4):
                        suml.append(A[l][i])
                    cur.append(sum(suml))
                sumP.append(np.array(cur))
            for k in Alist:
                plt.plot(x[0][0], sumP[list(Alist).index(k)], "--", label = f"{k}", lw = 0.7)


    if ferm == 1:
        plt.title(f"Neutron density distributions in Eu{Alist[0]} to Eu{Alist[-1]}")
    else:
        plt.title(f"Proton density distributions in Eu{Alist[0]} to Eu{Alist[-1]}")

    plt.ylabel(r"Density $\mathrm{fm}^{-3}$")
    plt.xlabel(r"x $\mathrm{fm}$")
    plt.legend(loc = "best")
    plt.xlim(0, 10)
    #plt.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data', f'dens_{Alist[0]}-{Alist[-1]}.jpg'), dpi = 500, bbox_inches = "tight", pad_inches = 0.06)
    plt.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data2',  "TB", f'dens_{Alist[0]}-{Alist[-1]}.jpg'), dpi = 500, bbox_inches = "tight", pad_inches = 0.06)
    print(f"Plot generated at '{os.path.join(os.getcwd(), 'unconstrained', 'data2', 'TB', f'dens_{Alist[0]}-{Alist[-1]}.jpg')}'")

else:
    dtheta = 2**(-10)
    
    for A in range(len(Alist)):
        if ferm == 1:
            den = nden[A]
        else:
            den = pden[A]
        densit = []
        pcord = []
        ccord = []

        if component == True:
            fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize=(15, 15))
            ax = ax.flatten()
            for axes in ax:
                axes.set(aspect='equal')
        else:
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.set(aspect='equal')

        azi = np.arange(0, 2*np.pi, dtheta)
        for angle in azi:
                for r in x[0][0]:
                    pcord.append((r, angle))
                    ccord.append((r*np.sin(angle), r*np.cos(angle)))

        densitinfo = []

        l = 0
        while l < 4:
            densitl = []
            densitlinfo = []
            for angle in azi:
                for r in x[0][0]:
                    fr = den[l][x[0][0].index(r)]
                    Ptheta = np.polynomial.legendre.legval(np.cos(angle), eval(f"({'0,'*l*2}1)"))
                    densitl.append(fr*Ptheta)
                    densitlinfo.append((fr*Ptheta, r, angle))
                    #print("{:5.0f}{:15.4E}{:10.4f}".format(l*2, r, angle*180/np.pi))

            densitinfo.append(densitlinfo)
            densit.append(densitl)
            xc, yc = [i[0] for i in ccord], [i[1] for i in ccord]

            if component == True:
                p = ax[l].tricontourf(xc, yc, densitl, 50, cmap = "plasma")
                cbar = fig.colorbar(p)
                cbar.ax.set_title(r"${fm}^{-3}$")
                ax[l].axis([-10, 10, -10, 10])
                ax[l].set_ylabel(r"z $\mathrm{fm}$")
                ax[l].set_xlabel(r"x $\mathrm{fm}$")
                ax[l].text(0.89, 0.94, r'$\lambda = ${:.2f}'.format(int(l*2)), horizontalalignment='center', fontsize = "10",
                        bbox=dict(facecolor='white', edgecolor='black', pad=3), transform=ax[l].transAxes)
            l += 1


        if component == True:
            #fig.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data', f'dens2dcomp_{Alist[A]}.jpg'), dpi = 500)
            fig.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data2', f'dens2dcomp_{Alist[A]}.jpg'), dpi = 500)
            print(f"Plot generated at '{os.path.join(os.getcwd(), 'unconstrained', 'data2', f'dens2dcomp_{Alist[A]}.jpg')}'")

        if component == False:
            tot = []
            for i in range(len(densit[0])):
                temp = []
                for l in densit:
                    temp.append(l[i])
                tot.append(sum(temp))
            p = ax.tricontourf(xc, yc, tot, 50, cmap = "plasma")

            ax.axis([-10, 10, -10, 10])
            ax.set_ylabel(r"z $\mathrm{fm}$")
            ax.set_xlabel(r"x $\mathrm{fm}$")
            cbar = fig.colorbar(p)
            cbar.ax.set_title(r"${fm}^{-3}$")
            if trace == 0 or trace == 2:
                #fig.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data', 'densconv_p', f'dens2d_{Alist[A]}.jpg'), dpi = 500)
                fig.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data2', 'densconv_p', f'dens2d_{Alist[A]}.jpg'), dpi = 500)
                print(f"Plot generated at '{os.path.join(os.getcwd(), 'unconstrained', 'data2', 'densconv_p', f'dens2d_{Alist[A]}.jpg')}'")


        if param == True:
            rhomax = []
            sumdensitinfo = []
            L = 0
            Q = 0
            
            for i in range(len(densitinfo[0])):
                sumdensitinfo.append((densitinfo[0][i][0]+densitinfo[1][i][0]+densitinfo[2][i][0]+
                                    densitinfo[3][i][0], densitinfo[0][i][1], densitinfo[0][i][2]))
                #print("Summing densities ({}/{})".format(int(i), int(len(densitinfo[0]))))

            for r in azi:
                curang = []
                for rho in sumdensitinfo:
                    #print("Searching theta ={:15.10f} ({}/{})".format(r*180/np.pi, 
                                                                    # int(sumdensitinfo.index(rho)+1), 
                                                                    # int(len(sumdensitinfo))))
                    if rho[2] == r:
                        curang.append(rho)
                maxind = [j[0] for j in curang].index(max([j[0] for j in curang]))
                rhomax.append(curang[maxind])
                dl = curang[maxind][1]*dtheta
                L += dl
                Q += curang[maxind][0]*dl
                
            rc, tc = [t[1] for t in rhomax], [u[2] for u in rhomax]
            xo = []
            zo = []
            for w in range(len(rc)):
                xo.append(rc[w]*np.sin(tc[w]))
                zo.append(rc[w]*np.cos(tc[w]))

            if component == False and (trace == 2 or trace == 1):
                ax.scatter(xo, zo, s = 0.3, marker = "x")
                #fig.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data', 'densconv_p', f'densmax_{Alist[A]}trace.jpg'), dpi = 500)
                fig.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data2', 'densconv_p', f'densmax_{Alist[A]}trace.jpg'), dpi = 500)
                print(f"Plot generated at '{os.path.join(os.getcwd(), 'unconstrained', 'data2','densconv_p', f'densmax_{Alist[A]}trace.jpg')}'")

            et = time.time()

            cden = sum([j[0] for j in den])
            #with open(os.path.join(os.getcwd(), 'unconstrained', 'data', 'densconv_p', 'bubble_param.txt'), "a") as convfile:
            with open(os.path.join(os.getcwd(), 'unconstrained', 'data2', 'densconv_p', 'bubble_param.txt'), "a") as convfile:
                convfile.write("{:<6.0f}{:<25.20f}{:<15.10f}{:<25.20f}{:<25.6f}\n".format(int(Alist[A]), Q/L, cden, 1-(cden/(Q/L)), (1-(cden/(Q/L)))*100))
            print(f"Data appended to '{os.path.join(os.getcwd(), 'unconstrained', 'data2', 'densconv_p', 'bubble_param.txt')}'")

            plt.close(fig = fig)

    if param == True:
        bp = []
        if ferm == 1:
            with open(os.path.join(os.getcwd(), 'unconstrained', 'data', 'densconv_n', 'bubble_param.txt'), "r") as convfile:
                convdata = convfile.readlines()
                for line in convdata[1:]:
                    bp.append(float(line[70:].strip()))
        else:
            #with open(os.path.join(os.getcwd(), 'unconstrained', 'data', 'densconv_p', 'bubble_param.txt'), "r") as convfile:
            with open(os.path.join(os.getcwd(), 'unconstrained', 'data2', 'densconv_p', 'bubble_param.txt'), "r") as convfile:
                convdata = convfile.readlines()
                for line in convdata[1:]:
                    bp.append(float(line[70:].strip()))

        fig2, ax2 = plt.subplots(figsize=(10, 10))
        ax2.plot(Alist, bp, "o-", color = "blue")
        ax2.set_xlabel("A")
        ax2.set_ylabel("DF")
        fermname = lambda x : "Neutron" if x == 1 else "Proton"
        ax2.set_title(f"Deformation factor at ground state for isotopes of Eu ({fermname(ferm)})")
        Alistlabel=[]
        for A in Alist:
            Alistlabel.append(f"{A}")
        ax2.set_xticks(Alist, Alistlabel)
        #fig2.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data', f'df.jpg'), dpi = 500)
        fig2.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data2', f'df.jpg'), dpi = 500)
