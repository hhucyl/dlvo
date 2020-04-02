import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Circle

fn = "test_swi_ini.txt"
f = open(fn,"w")
nu = 0.01
Re = 1e2

f.write(str(nu)+"\tnu\n")
f.write(str(Re)+"\tRe\n")

R = 10.0
RR = 50.0
tpl = 3
rwl = 6*R
fpl = 6*R+tpl
bbl = 2*R+2
d = 2.0
vmax = Re*nu/(rwl+fpl)*1.5
print("vmax ",vmax)
f.write(str(vmax)+"\tmax vel\n")
f.write(str(R)+"\tR\n")
f.write(str(RR)+"\tRR\n")
f.write(str(d)+"\td hydraulic path\n")
f.write(str(rwl)+"\trwl the space in the random walk layer\n")
f.write(str(fpl)+"\tfpl the space in the fine particle layer\n")
f.write(str(bbl)+"\tbbl the space in the bottom buffer layer\n")
f.write(str(tpl)+"\ttpl the space in the top buffer layer\n")

Pnx = 3
Pny = 3
fpn = 20
pdx = 60
pdy = 60

f.write(str(Pnx)+"\tPnx big Particle num in x\n")
f.write(str(Pny)+"\tPny big Particle num in y\n")
f.write(str(fpn)+"\tSmall Particle num\n")
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
Gn = 0.0
Kt = 0.0
Mu = 0.0
Eta = 0
Beta = 0
# A = 3e-20
# kappa = 1e9
# fi = 15*5.4
# Z = 9.22e-11*np.tanh(fi/103.0)*np.tanh(fi/103.0)
A = 0
kappa = 0
fi = 0
Z = 0
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
GY = 2 #GY= 1 has bug

f.write(str(GX)+"\tGX\n")
f.write(str(GY)+"\tGY\n")

A = A/((ratiol/ratiot)*(ratiol/ratiot)*ratiom)
kappa = kappa*ratiol
Z = Z*ratiot*ratiot/(ratiol*ratiom)
print("A ",A)
print("kappa ",kappa)
print("Z ",Z)


nx = np.ceil(Pnx*(RR*2 + pdx))

H =  rwl + fpl + 1
sy =  (Pny-1)*(np.sqrt(3)*RR+pdy) + 2*RR + bbl
ny = np.ceil(sy+H)
f.write(str(int(nx))+"\tnx\n")
f.write(str(int(ny))+"\tny\n")
print("sy ",sy)
pl = 1.0*nx/(fpn*1.0)
print(pl)
print("nx ", nx)
print("ny ", ny)
ll = 1.2*R
f.write(str(ll)+"\tll buffer for add\n")

ppx = np.zeros(fpn+Pnx*Pny)
ppy = np.zeros(fpn+Pnx*Pny)
ptag = np.zeros(fpn+Pnx*Pny)
pr = np.zeros(fpn+Pnx*Pny)
fig = plt.figure()
ax = fig.add_subplot(111)
pn = 0

# for i in range(fpn):
#     px = 0.5*pl + i*pl
#     py = np.random.uniform(sy + rwl + ll, ny-1 -R - ll)
#     ppx[pn] = px
#     ppy[pn] = py
#     pr[pn] = R
#     ptag[pn] = 1 
#     pn = pn+1
#     cir = Circle(xy=(px,py),radius=R,alpha=0.2)
#     ax.add_patch(cir)

for i in range(fpn):
    
    while True:
        rand = np.random.rand(2)
        px = np.random.uniform(ll, nx-1-ll)
        py = np.random.uniform(sy + rwl + ll, ny-1 -ll- tpl)
        nn = 0
        overlap = False
        for ii in range(i):
            dd = np.sqrt((px-ppx[ii])**2+(py-ppy[ii])**2)
            if dd<3*R:
                overlap = True
                break
        if not overlap:
            break

    ppx[i] = px
    ppy[i] = py
    pr[pn] = R
    ptag[pn] = 1 
    pn = pn+1

    cir = Circle(xy=(ppx[i],ppy[i]),radius=R,alpha=0.2)
    ax.add_patch(cir)
print("sy ",sy)
for j in range(Pny):
    py = j*(np.sqrt(3)*RR+pdy) + RR + 1 + bbl
    temp = Pnx if j%2==0 else Pnx+1
    temp = Pnx
    for i in range(temp):
        px = 0.5*pdx + RR + i*(2*RR+pdx) if j%2==0 else i*(2*RR+pdx)

        ppx[pn] = px
        if j==Pny-1 and i == 1:
            ppy[pn] = py + d
        else:
            ppy[pn] = py
        pr[pn] = RR
        ptag[pn] = -1 
        pn = pn+1

        cir = Circle(xy=(px,py),radius=RR,alpha=0.2)
        ax.add_patch(cir)

sy = np.max(ppy[np.where(ptag == -1)]) + RR


plt.plot([0,nx-1],[0,0],'k')
plt.plot([0,nx-1],[ny-1,ny-1],'k')
plt.plot([0,0],[0,ny-1],'k')
plt.plot([nx-1,nx-1],[0,ny-1],'k')

plt.plot([0,nx-1],[sy,sy],'r',label='sy')
plt.plot([0,nx-1],[bbl,bbl],'b',label='bbl')
plt.plot([0,nx-1],[rwl+sy,rwl+sy],'b',label='rwl')
plt.title(fpn)
plt.axis('equal')
plt.savefig("fine.png")
plt.clf()
cutoff = 1.57e-10/ratiol
d = np.arange(cutoff,1e-2,ratiol)
Fvdw = A/(6*d**2)*R*RR/(R+RR)
Fe = kappa*R*RR/(R+RR)*Z*np.exp(-kappa*d)
plt.plot(d,Fvdw,'r',label='Vdw Force')
plt.plot(d,Fe,'b',label='Elec Force')
plt.plot(d,Fvdw-Fe,'k',label='Total Force')
plt.legend()
plt.savefig("Force.png")

f.write(str(sy)+'\tsy\n')
f.write(str(fpn+Pnx*Pny)+'\tpinit\n')
for i in range(fpn+Pnx*Pny):
    f.write(str(ppx[i])+' '+str(ppy[i])+' '+str(pr[i])+' '+str(ptag[i])+'\n')
f.close()



