import numpy as np
import matplotlib.pyplot as plt

RR = 100
R = 10
Pnx = 5
Pny = 5
pnx = 10
pny = 10
dx = 0
dy = 0
ppl = 2*RR
Nx = Pnx*(dx+2*RR)
pl = Nx/pnx;
Ny = (Pny-1)*(3**0.5*RR+dy)+2*RR+ppl+(pny-1)*pl+6*R
sy = (Pny-1)*(3**0.5*RR+dy)+2*RR+ppl

px = 0
py = 0

fig = plt.gcf()
ax = fig.gca()
for j in range(Pny):
	py = 0 + j*(3**0.5*RR+dy) + RR
	temp = Pnx if np.mod(j,2)==0 else Pnx+1
	for i in range(temp):
		if np.mod(j,2)==0:
			px = 0.5*dx+RR + i*(2*RR+dx)
		else:
			px = i*(2*RR+dx)
		circle = plt.Circle((px,py),0.8*RR,fill=False)
		ax.add_patch(circle)

for j in range(pny):
	py = sy + j*pl
	for i in range(pnx):
		px = 0.5*pl + i*pl
		circle = plt.Circle((px,py),R,fill=False)
		ax.add_patch(circle)
ax.set_aspect('equal', adjustable='datalim')
ax.plot([0,0],[0, Ny],'k')
ax.plot([Nx,Nx],[0, Ny],'k')
ax.plot([0,Nx],[0, 0],'k')
ax.plot([0,Nx],[Ny, Ny],'k')

plt.show()
