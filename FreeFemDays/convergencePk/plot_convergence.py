import numpy as np 
import re
from collections import defaultdict
import matplotlib.pyplot as plt 

filename = "Poisson_2DC.txt"

# dictionnaire pour stocker
data = defaultdict(list)

current_p = None# dictionnaire pour stocker
data = defaultdict(list)

current_p = None

with open(filename, "r") as f:
    for line in f:
        line = line.strip()
        
        # si la ligne contient "P=..."
        match_p = re.match(r"P=(\d+)", line)
        if match_p:
            current_p = int(match_p.group(1))
            continue

        # si la ligne contient "Npoint(...)=... Error(relative)=..."
        match_data = re.match(r"Npoint\(npplo\)=(\d+)\s+Error\(relative\)=\s*([0-9.eE+-]+)", line)
        if match_data and current_p is not None:
            npoint = int(match_data.group(1))
            error = float(match_data.group(2))
            data[current_p].append((npoint, error))

for p, values in data.items():
    print(f"P={p}")
    for npoint, error in values:
        print(f"  N={npoint}, Error={error}")
plt.figure()
colors=["blue","orange","green","red","black","brown","magenta"]
cst=[25,30,16.85,8,8.33,3,2.5]
for p in data.keys():
    npoints = np.array([n for n, err in data[p]])
    errors = np.array([err for n, err in data[p]])

    print("P=",p)
    print("npoints=",npoints)
    print("errors=",errors)
    print("Pente=",(np.log(errors[1:])-np.log(errors[:-1]))/(np.log(npoints[1:])-np.log(npoints[:-1])) )
    plt.loglog(npoints,errors,label="P"+str(p),marker="x",color=colors[p-1])
    
    #plt.loglog(npoints,cst[p-1]*(1/npoints)**(p+1),label = rf"$1/h^{{{p+1}}}$",linestyle="dashed",color=colors[p-1])

plt.legend(ncol=2,fontsize=9)
#plt.yscale("log")
C1 = 1e-6 
C2 = 1.0
C3 = 10.0
C4 = 100.0

'''plt.loglog(npoints,25*(1/npoints)**2,label="1/h^2",linestyle="dashed",color="blue")
plt.loglog(npoints,0.1*(1/npoints)**3,label="1/h^3",linestyle="dashed",color="orange")
plt.loglog(npoints,0.001*(1/npoints)**4,label="1/h^4",linestyle="dashed",color="green")
plt.loglog(npoints,80*(1/npoints)**5,label="1/h^5",linestyle="dashed",color="red")
plt.loglog(npoints,80*(1/npoints)**6,label="1/h^6",linestyle="dashed",color="black")
plt.loglog(npoints,80*(1/npoints)**7,label="1/h^7",linestyle="dashed",color="brown")'''

plt.title("Relative Error")




plt.xlabel("Nbr of points")
plt.ylabel("Relative Error")
plt.savefig("Poisson_2DC_sanspente.pdf",format="pdf")
plt.show()