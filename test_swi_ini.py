import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Circle

fn = "test_swi_ini.txt"
f = open(fn,"w")
nu = 0.01
Re = 1e3

f.write(str(nu)+"\tnu\n")
f.write(str(Re)+"\tRe\n")

R = 10.0
RR = 50.0
ppl = 2*R
pppl = 2*R

f.write(str(R)+"\tR\n")
f.write(str(RR)+"\tRR\n")
f.write(str(ppl)+"\tppl the space in the bottom\n")
f.write(str(pppl)+"\tpppl the space in the up\n")

Pnx = 3
Pny = 3
pnx = 10
pny = 8
pdx = 60
pdy = 60

f.write(str(Pnx)+"\tPnx big Particle num in x\n")
f.write(str(Pny)+"\tPny big Particle num in y\n")
f.write(str(pnx)+"\tpnx small Particle num in x\n")
f.write(str(pny)+"\tpny small Particle num in y\n")
f.write(str(pdx)+"\tpdx gap in x\n")
f.write(str(pdy)+"\tpdy gap in y\n")

ratiol = 1e-6
ratiot = 1e-8
ratiom = 1e-15

f.write(str(ratiol)+"\tratiol\n")
f.write(str(ratiot)+"\tratiot\n")
f.write(str(ratiom)+"\tratiom\n")

rho = 1.0
rhos = 2.7
Ga = 100.0

f.write(str(rho)+"\trho\n")
f.write(str(rhos)+"\trhos\n")
f.write(str(Ga)+"\tGa\n")

Kn = 50
Gn = 0
Kt = 0
Mu = 0
Eta = 0
Beta = 0
A = 3e-20
kappa = 1e9
Z = 3.97e-11
D = 2

f.write(str(Kn)+"\tKn\n")
f.write(str(Gn)+"\tGn\n")
f.write(str(Kt)+"\tKt\n")
f.write(str(Mu)+"\tMu\n")
f.write(str(Eta)+"\tEta\n")
f.write(str(Beta)+"\tBeta\n")
f.write(str(A)+"\tA\n")
f.write(str(kappa)+"\tkappa\n")
f.write(str(Z)+"\tZ\n")
f.write(str(D)+"\tD\n")

GX = 3
GY = 2

f.write(str(GX)+"\tGX\n")
f.write(str(GY)+"\tGY\n")

A = A/((ratiol/ratiot)*(ratiol/ratiot)*ratiom)
kappa = kappa*ratiol
Z = Z*ratiot*ratiot/(ratiol*ratiom)
print("A ",A)
print("kappa ",kappa)
print("Z ",Z)


nx = np.ceil(Pnx*(RR*2 + pdx))
ppl = 2*R
pl = nx/pnx
H =  ppl + (pny-1)*pl + pppl + 1
ll = R*1.2



ppx = np.zeros(pnx*pny)
ppy = np.zeros(pnx*pny)

fig = plt.figure()
ax = fig.add_subplot(111)
N = 0
for i in range(pnx*pny):
    
    while True:
        rand = np.random.rand(2)
        px = ll+rand[0]*(nx-1-ll-ll)
        py = ll+ppl+rand[1]*(H-pppl-ll-ll-ppl)
        nn = 0
        overlap = False
        for ii in range(i):
            d = np.sqrt((px-ppx[ii])**2+(py-ppy[ii])**2)
            if d<3*R:
                overlap = True
                break
        if not overlap:
            break

    ppx[i] = px
    ppy[i] = py
    if(ppy[i]>0.5*H):
        N = N+1
    cir = Circle(xy=(ppx[i],ppy[i]),radius=R,alpha=0.2)
    ax.add_patch(cir)



plt.plot([0,nx-1],[0,0],'k')
plt.plot([0,nx-1],[H,H],'k')
plt.plot([0,0],[0,H],'k')
plt.plot([nx-1,nx-1],[0,H],'k')
plt.plot([0,nx-1],[0.5*H,0.5*H],'b')
plt.title(N)
plt.axis('equal')
plt.savefig("fine.png")
plt.clf()
cutoff = 1.57e-10/ratiol
d = np.arange(cutoff,1e-2,ratiol)
Fvdw = A/(6*d**2)*R*RR/(R+RR)
Fe = kappa*R*RR/(R+RR)*Z*np.exp(-kappa*d)
plt.plot(d,Fvdw,'r')
plt.plot(d,Fe,'b')
plt.plot(d,Fvdw-Fe,'k')
plt.savefig("Force.png")

f.write(str(N)+'\tpinit\n')
for i in range(pnx*pny):
    f.write(str(ppx[i])+' '+str(ppy[i])+'\n')
f.close()



