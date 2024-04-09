import numpy as np
import matplotlib.pyplot as plt
import os 
import random
import matplotlib.pylab as pl
import matplotlib.patches as patches

ferm = 2 # fermion 1=neutron 2=proton
par = [1, 1] #parity  [0=both, 1=single, 1=+ 2=-]
lock = False
maxE = 2
minE = -90
Alist = [179, 181, 183, 185, 187, 189]
select = False
orbitals = ["3s 1/2"]

levelList = {}
palette = []
edgePoints = []   
fermiPoints = []
cmap = pl.cm.hsv(np.linspace(0.0,1.0,30))

for A in Alist:
    Afile = os.path.join(os.getcwd(), "unconstrained", "data", f"{A}")
    searchFile = os.path.join(Afile, f"{A}u_bulk.txt")
    with open(searchFile, "r") as searchData:
        data = searchData.readlines()
        target = f" Etot + Ecm                                            {data[-1][:12]}\n"
        for line in data:
            if line == target:
                minB = round(float(data[data.index(line)-5][7:].strip()), 2)
                if minB == 0:
                    minB = 0
    
    dataFile = os.path.join(os.getcwd(), "groundstates.txt")
    fermi = []
    with open(dataFile, "r") as ground:
        data = ground.readlines()
        for line in data:
            if line[0:4].strip() == str(A):
                fermi = [float(line[105:114].strip()), float(line[115:125].strip())]

    plt.plot([A-0.5, A+0.5], [fermi[ferm-1], fermi[ferm-1]], "-", lw = 0.75, marker = "", c = 'black')
    if list(Alist).index(A) == 0:
        fermiPoints.append([A+0.5, fermi[ferm-1]])
        if ferm == 1:
            lamb = r"$\lambda_N$"
        else:
            lamb = r"$\lambda_P$"
        plt.annotate(lamb, xy=(A-0.5, fermi[ferm-1]), xytext =(A-0.75, fermi[ferm-1]), fontsize = 4, c = "black", 
                     verticalalignment = "center", horizontalalignment = "center")
    elif list(Alist).index(A) == -1:
        fermiPoints.append([A-0.5, fermi[ferm-1]])
    else:
        fermiPoints.append([A-0.5, fermi[ferm-1]])
        fermiPoints.append([A+0.5, fermi[ferm-1]])

    splFile = os.path.join(Afile, f"initb.{minB}", "spl.out")

    title = "it ip   m    i      ecan       eqp      dcan        vcan\n"
    levels = {}
    with open(splFile, "r") as splData:
        data = splData.readlines()
        ind = list(np.where(np.array(data) == title)[0])
        a, b = 0, 1
        no = len(ind)
        iter = 1
        if par[0] == 1:
            while iter <= no:
                start = ind[a]
                if float(data[start+1][1].split()[0]) == ferm and float(data[start+1][3].split()[0]) == par[1]:
                    if iter == no:
                        end = ""
                    else:
                        end = ind[b]
                    E = float(data[start+1][16:26].split()[0])
                    levels[E] = eval(f"data[start:{end}]") 
                else:
                    pass
                iter += 1
                a += 1
                b += 1
        else:
            while iter <= no:
                start = ind[a]
                if float(data[start+1][1].split()[0]) == ferm:
                    if iter == no:
                        end = ""
                    else:
                        end = ind[b]
                    E = float(data[start+1][16:26].split()[0])
                    levels[E] = eval(f"data[start:{end}]") 
                else:
                    pass
                iter += 1
                a += 1
                b += 1
            
    sortE = sorted(levels.keys(), reverse = True)

    label = lambda x : "_neu" if x == 1 else "_pro"

    if os.path.exists(os.path.join(Afile, f"spl{label(ferm)}full.out")) == True:
        print(f"Sorted version already exists at '{os.path.join(Afile, f'spl{label(ferm)}full.out')}'\n")
    else:
        with open(os.path.join(Afile, f"spl{label(ferm)}full.out"), "a") as splWrite:
            for e in list(sortE):
                for line in list(levels.values())[list(levels.keys()).index(e)]:
                    splWrite.write(line)
        print(f"'{splFile}' sorted as '{os.path.join(Afile, f'spl{label(ferm)}full.out')}'\n")
    
    Ecrop = []
    levelName = []
    amp = []
    ampProb = []
    vcan = []

    eA = np.array(sortE, dtype = float)
    A1 = maxE > eA
    A2 = eA > minE
    match = A1*1 + A2*1
    for e in list(eA[match == 2]):
        if select == False:
            Ecrop.append(e)
            v = float(list(levels.values())[list(levels.keys()).index(e)][1][47:].strip())
            vcan.append(v)
            id = str(list(levels.values())[list(levels.keys()).index(e)][1][3:10].strip())
            levelTemp = []
            ampTemp = []
            for level in list(levels.values())[list(levels.keys()).index(e)][3:]:
                levelTemp.append(level[7:13].strip())
                ampTemp.append(float(level[34:47].strip()))
            levelName.append([levelTemp[ampTemp.index(max(ampTemp))], id])
            amp.append(max(ampTemp))
            ampProb.append(max(ampTemp)*v)
        else:
            v = float(list(levels.values())[list(levels.keys()).index(e)][1][47:].strip())
            id = str(list(levels.values())[list(levels.keys()).index(e)][1][3:10].strip()) 
            for level in list(levels.values())[list(levels.keys()).index(e)][3:]:
                if level[7:13].strip() in orbitals:
                    Ecrop.append(e)
                    vcan.append(v)
                    levelName.append([level[7:13].strip(), id])
                    amp.append(float(level[34:47].strip()))
                    ampProb.append(float(level[34:47].strip())*v)
                else:
                    pass
    
    if os.path.exists(os.path.join(Afile, f"spl{label(ferm)}.out")) == False:
        with open(os.path.join(Afile, f"spl{label(ferm)}.out"), "a") as splWrite:
            for e in list(eA[match == 2]):
                Ecrop.append(e)
                for line in list(levels.values())[list(levels.keys()).index(e)]:
                    splWrite.write(line)           

        with open(os.path.join(Afile, f"spl{label(ferm)}.out"), "a") as splWrite:
            splWrite.write(f"{Ecrop}\n")
            splWrite.write(f"{levelName}\n")
            splWrite.write(f"{amp}\n")
            splWrite.write(f"{ampProb}\n")

        print(f"Cropped version of '{splFile}' created at '{os.path.join(Afile, f'spl{label(ferm)}.out')}'\n")
    else:
        print(f"Cropped version of '{splFile}' already exists at '{os.path.join(Afile, f'spl{label(ferm)}.out')}'\n")

    for L in levelName:
        if L[0] not in levelList.keys():
            levelList[L[0]] = L[1]
            repeat = True
            while repeat == True:
                c = random.randint(0,29)
                if c not in palette:
                    palette.append(c)
                    repeat = False
        else:
            pass

    set = []
    count = 0
    for l in levelName:
        colour = cmap[palette[list(levelList.keys()).index(l[0])]]
        if select == True:
            high = amp.index(max(amp))
            if count == high:
                plt.plot([A-0.5, A+(0.5*vcan[count])], [Ecrop[count], Ecrop[count]], 
                 "-", lw = 0.75, marker = "|", markersize = 2.5, c = colour)
                plt.plot([A+(0.5*vcan[count]), A+0.5], [Ecrop[count], Ecrop[count]], 
                 "-", lw = 0.3, marker = "|", markersize = 2.5, c = "grey")
                plt.annotate("{:.2f}".format(vcan[count]), xy=(A, Ecrop[count]), 
                        xytext =(A, Ecrop[count]+0.2), 
                        c = colour, horizontalalignment='center', 
                        verticalalignment = 'center',  fontsize = 5)
                
                if list(Alist).index(A) == 0:
                    set.append((0, [A+0.5, Ecrop[count]], (l)))
                elif list(Alist).index(A) == len(list(Alist))-1:    
                    set.append(([A-0.5, Ecrop[count]], 0, (l)))
                else:
                    set.append(([A-0.5, Ecrop[count]], [A+0.5, Ecrop[count]], (l)))
            else:
                pass
            
            count += 1 
        
        else:
            plt.plot([A-0.5, A+(0.5*vcan[levelName.index(l)])], [Ecrop[levelName.index(l)], Ecrop[levelName.index(l)]], 
                 "-", lw = 0.75, marker = "|", markersize = 2.5, c = colour)
            plt.plot([A+(0.5*vcan[levelName.index(l)]), A+0.5], [Ecrop[levelName.index(l)], Ecrop[levelName.index(l)]], 
                 "-", lw = 0.3, marker = "|", markersize = 2.5, c = "grey")
        
            if list(Alist).index(A) == 0:
                set.append((0, [A+0.5, Ecrop[levelName.index(l)]], (l)))
            elif list(Alist).index(A) == len(list(Alist))-1:    
                set.append(([A-0.5, Ecrop[levelName.index(l)]], 0, (l)))
            else:
                set.append(([A-0.5, Ecrop[levelName.index(l)]], [A+0.5, Ecrop[levelName.index(l)]], (l)))
          
    edgePoints.append(set)
        

"""     for p in vcan:
        if p == 0:
            pass
        else:
            valence = vcan.index(p)
            break
    plt.plot(A, Ecrop[valence], marker = "D", c = cmap[palette[list(levelList.keys()).index(levelName[valence][0])]], 
             markersize = 3) """

labelled = []
for iso in edgePoints:
    for level in iso:
        name = level[2][0]
        if name not in labelled:
            if edgePoints.index(iso) == len(edgePoints)-1:
                plt.annotate(f"{name}", xy=(Alist[edgePoints.index(iso)]-0.5, level[0][1]), 
                    xytext =(Alist[edgePoints.index(iso)]-0.75, level[0][1]), 
                    c = cmap[palette[list(levelList.keys()).index(name)]], horizontalalignment='center',
                    verticalalignment = 'center',  fontsize = 4)
            else:
                plt.annotate(f"{name}", xy=(Alist[edgePoints.index(iso)]-0.5, level[1][1]), 
                    xytext =(Alist[edgePoints.index(iso)]-0.75, level[1][1]), 
                    c = cmap[palette[list(levelList.keys()).index(name)]], horizontalalignment='center',
                    verticalalignment = 'center',  fontsize = 4)
            labelled.append(name)
        else:
            pass

s = 1
while s < len(edgePoints):
    nextIso = edgePoints[s]
    curIso = edgePoints[s-1]
    for nextlev in nextIso:
        iden = nextlev[2]
        P2 = nextlev[0]
        for curlev in curIso:
            if curlev[2] == iden:
                P1 = curlev[1]
                plt.plot([P1[0], P2[0]], [P1[1], P2[1]], "--", marker = "", lw = 0.5, 
                         c = cmap[palette[list(levelList.keys()).index(iden[0])]])
            else:
                pass

    if s == len(edgePoints) - 1:
        for nextlev in edgePoints[s]:
            iden = nextlev[2]
            parity = lambda x : "+" if x == 1 else "-"
            plt.annotate(f"{int(int(iden[1][-2:].strip())*2-1)}/2{parity(int(iden[1][0].strip()))}", 
                        xy=(Alist[s]+0.5, nextlev[0][1]), xytext =(Alist[s]+0.75, nextlev[0][1]), 
                        c = cmap[palette[list(levelList.keys()).index(iden[0])]], horizontalalignment='center', 
                        verticalalignment = 'center', fontsize = 4)
    else:
        pass
    s += 1

p = 0
while p < len(fermiPoints)-2:
    plt.plot([fermiPoints[p][0], fermiPoints[p+1][0]], [fermiPoints[p][1], fermiPoints[p+1][1]], '--',lw = 0.5, c = 'black')
    p += 2

legend = []
for n in levelList.keys():
    legend.append(patches.Patch(color=cmap[palette[list(levelList.keys()).index(n)]], label=f"{n}"))
plt.legend(handles = legend, loc='center right', bbox_to_anchor=(1.2, 0.5), ncol=1, fontsize = 8)
if maxE > 0:
    plt.axhline(y = 0, color = 'silver', linestyle = '-', lw = 0.4)
plt.xlabel(r"A(N)")
Alistlabel=[]
for A in Alist:
    Alistlabel.append(f"{A}({A-63})")
plt.ylabel(r"$E_{\mathrm{can}}$ ($\mathrm{MeV}$)")
plt.xticks(Alist, Alistlabel)
fermname = lambda x : "neutron" if x == 1 else "proton"
parname = lambda x : "+ve" if x == 1 else "-ve"
if lock == True:
    plt.ylim(minE, maxE)
if par[0] == 1:
    plt.title(f"Single {fermname(ferm)} energies of Eu isotopes ({parname(par[1])} parity)")
else:
    plt.title(f"Single {fermname(ferm)} energies of Eu isotopes (All parities)")
plt.tick_params(axis = "y", right = True)
if select == False:
    plt.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data', f'spl_{Alist[0]}-{Alist[-1]}.jpg'), dpi = 500, bbox_inches = "tight", pad_inches = 0.06)
    print(f"Plot generated at '{os.path.join(os.getcwd(), 'unconstrained', 'data', f'spl_{Alist[0]}-{Alist[-1]}.jpg')}'")
else:
    plt.savefig(os.path.join(os.getcwd(), 'unconstrained', 'data', f'spl_{Alist[0]}-{Alist[-1]}select.jpg'), dpi = 500, bbox_inches = "tight", pad_inches = 0.06)
    print(f"Plot generated at '{os.path.join(os.getcwd(), 'unconstrained', 'data', f'spl_{Alist[0]}-{Alist[-1]}select.jpg')}'")