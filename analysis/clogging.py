import numpy as np
import matplotlib.pyplot as plt

RR = 100
RRR = 3*RR
Pnx = 50
Pny = 49
dx = 0
dy = 0
ppl = 2*RR
Nx = Pnx*(dx+2*RR)
Ny = (Pny-1)*(3**0.5*RR+dy)+2*RR

px = 0
py = 0

fig = plt.gcf()
ax = plt.subplot(121)
for j in range(Pny):
	py = 0 + j*(3**0.5*RR+dy) + RR
	temp = Pnx if np.mod(j,2)==0 else Pnx+1
	for i in range(temp):
		if np.mod(j,2)==0:
			px = 0.5*dx+RR + i*(2*RR+dx)
		else:
			px = i*(2*RR+dx)
		circle = plt.Circle((px,py),RR,fill=False)
		ax.add_patch(circle)

ax.set_aspect('equal', adjustable='datalim')
ax.plot([0,0],[0, Ny],'k')
ax.plot([Nx,Nx],[0, Ny],'k')
ax.plot([0,Nx],[0, 0],'k')
ax.plot([0,Nx],[Ny, Ny],'k')

fig = plt.gcf()
ax = plt.subplot(122)
ppx = 0.5*Nx
ppy = 0.5*Ny
for j in range(Pny):
	py = 0 + j*(3**0.5*RR+dy) + RR
	temp = Pnx if np.mod(j,2)==0 else Pnx+1
	for i in range(temp):
		if np.mod(j,2)==0:
			px = 0.5*dx+RR + i*(2*RR+dx)
		else:
			px = i*(2*RR+dx)
		if (px-ppx)**2+(py-ppy)**2<RRR*RRR:
			circle = plt.Circle((px,py),RR,fill=False)
			ax.add_patch(circle)

ax.set_aspect('equal', adjustable='datalim')
ax.plot([0,0],[0, Ny],'k')
ax.plot([Nx,Nx],[0, Ny],'k')
ax.plot([0,Nx],[0, 0],'k')
ax.plot([0,Nx],[Ny, Ny],'k')


plt.show()
