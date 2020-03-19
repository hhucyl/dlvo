import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Circle

fn = "test_w_ini.txt"
f = open(fn,"w")
nu = 0.01
Re = 1e1

f.write(str(nu)+"\tnu\n")
f.write(str(Re)+"\tRe\n")

R = 10.0
RR = 12.0
d = 2.0
ppl = 2*R
pppl = 2*R

f.write(str(R)+"\tR\n")
f.write(str(RR)+"\tRR\n")
f.write(str(d)+"\td hydraulic path\n")
f.write(str(ppl)+"\tppl the space in the bottom\n")
f.write(str(pppl)+"\tpppl the space in the top\n")

Pnx = 3
Pny = 3
pnx = 0
pny = 0
mag = 0.01
pdx = 2*mag*RR
pdy = 2*mag*RR

f.write(str(Pnx)+"\tPnx big Particle num in x\n")
f.write(str(Pny)+"\tPny big Particle num in y\n")
f.write(str(pnx)+"\tpnx small Particle num in x\n")
f.write(str(pny)+"\tpny small Particle num in y\n")
f.write(str(pdx)+"\tpdx gap in x\n")
f.write(str(pdy)+"\tpdy gap in y\n")
f.write(str(mag)+"\trandom ratio in fix top\n")

nx = int(np.ceil(Pnx*(RR*2 + pdx)))
pl = 1000000 if pnx==0 else nx/pnx

# H =  ppl + (pny-1)*pl + pppl + 1
H =  50+1
sy = Pny*(pdy+2*(RR-d)) + 1 +R
ny = int(np.ceil(sy+H))
f.write(str(nx)+"\tnx\n")
f.write(str(ny)+"\tny\n")
f.write(str(H)+"\tH\n")
f.write(str(pl)+"\tpl\n")

ratiol = 1e-6
ratiot = 1e-8
ratiom = 1e-15

f.write(str(ratiol)+"\tratiol\n")
f.write(str(ratiot)+"\tratiot\n")
f.write(str(ratiom)+"\tratiom\n")

rho = 1.0
rhos = 2.7
Ga = 10.0

f.write(str(rho)+"\trho\n")
f.write(str(rhos)+"\trhos\n")
f.write(str(Ga)+"\tGa\n")

# Kn = 50
# Gn = -0.3
# Kt = 10
# Mu = 0.4
# Eta = 0
# Beta = 0
# A = 3e-20
# kappa = 1e9
# fi = 15*5.4
# Z = 9.22e-11*np.tanh(fi/103.0)*np.tanh(fi/103.0)
# D = 2

Kn = 50
Gn = 0
Kt = 0
Mu = 0
Eta = 0
Beta = 0
A = 0
kappa = 0
fi = 15*5.4
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
GY = 2

f.write(str(GX)+"\tGX\n")
f.write(str(GY)+"\tGY\n")


A = A/((ratiol/ratiot)*(ratiol/ratiot)*ratiom)
kappa = kappa*ratiol
Z = Z*ratiot*ratiot/(ratiol*ratiom)
print("A ",A)
print("kappa ",kappa)
print("Z ",Z)



print("H", H)
ll = R*1.2



ppx = np.zeros(Pnx*Pny+pnx*pny)
ppy = np.zeros(Pnx*Pny+pnx*pny)
ptag = np.zeros(Pnx*Pny+pnx*pny)
pr = np.zeros(Pnx*Pny+pnx*pny)

fig = plt.figure()
ax = fig.add_subplot(111)
nn = 0

for j in range(Pny):
    py0 = 0.0 + j*(2*RR+pdy) + RR + 1 + pdy
    for i in range(Pnx):
        px = i*(2*RR+pdx) + pdx/2.0 + RR
        r = np.random.uniform(-mag, mag)
        px = px + r*RR
        py = py0 + r*RR
        
        ppx[nn] = px
        ppy[nn] = py
        pr[nn]  = RR
        ptag[nn] = 1#-1 fix 1 move
        nn = nn + 1
        cir = Circle(xy=(px,py),radius=RR-d,alpha=0.2)
        cir1 = Circle(xy=(px,py),color='r',radius=RR,alpha=0.1)
        ax.add_patch(cir) 
        ax.add_patch(cir1) 
#move ignore here
N = 0
for i in range(pnx*pny):
    print 1
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
    ppy[i] = py + sy
    pr[i]  = R
    ptag[i] = 1    

    if(ppy[i]>0.5*H):
        N = N+1
    cir = Circle(xy=(ppx[i],ppy[i]),radius=R-d,alpha=0.2)
    ax.add_patch(cir)

print('N ',N)

plt.plot([0,nx-1],[0,0],'k')
plt.plot([0,nx-1],[ny-1,ny-1],'k')
plt.plot([0,0],[0,ny-1],'k')
plt.plot([nx-1,nx-1],[0,ny-1],'k')
plt.plot([0,nx-1],[0.5*H+sy,0.5*H+sy],'b')
plt.plot([0,nx-1],[sy,sy],'b')
plt.title(N)
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

f.write(str(N)+'\tpinit\n')
f.write(str(sy)+'\tsy\n')

f.write(str(Pnx*Pny+pnx*pny)+'\tTotal Num\n')

for i in range(pnx*pny+Pnx*Pny):
    f.write(str(ppx[i])+' '+str(ppy[i])+' '+str(pr[i])+' '+str(ptag[i])+'\tpx, py, pr, ptag\n')
f.close()



