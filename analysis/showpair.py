import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import sys

prefix = "/media/pzhang/Ubuntu 16.0/"
prefix = "/home/user/dlvo/"
prefix = prefix + "test_swi_"
num = 20
R = 10
ratio = 8
RR = ratio*R
ratiol = 10e-6/(2.0*R)
ratiot = 1/1   
A = 3e-20/((ratiol/ratiot)*(ratiol/ratiot))

name = prefix + str(num).zfill(4) + ".h5"
f = h5.File(name)
Np = f['Np'][0]
Nx = f['Nx'][0]
Ny = f['Ny'][0]

pos = np.array(f['Pposition'])
pfc = np.array(f['PForce'])
pfh = np.array(f['PForceh'])
isf = np.array(f['PIsFree'])
posx = pos[0:3*Np-2:3]
posy = pos[1:3*Np-1:3]
pfcx = pfc[0:3*Np-2:3]
pfcy = pfc[1:3*Np-1:3]
pfhx = pfh[0:3*Np-2:3]
pfhy = pfh[1:3*Np-1:3]
fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(Np):
	if(isf[i]<0):
		cir = Circle(xy=(posx[i],posy[i]),radius=RR,alpha=0.2)
	else:
		cir = Circle(xy=(posx[i],posy[i]),radius=R,alpha=0.2)
	ax.add_patch(cir)
	ax.text(posx[i],posy[i],str(i))
ax.quiver(posx,posy,pfcx,pfcy,scale=0.2)
d = R+RR - np.sqrt((posx[6]-posx[19])**2+(posy[6]-posy[19])**2)
vd = 0.0005
Fvdw = A/(6*vd*vd)*R*RR/(R+RR)
print(Fvdw)
print(d,5*d)
print(pfcx[6],pfcy[6],np.sqrt(pfcx[6]**2+pfcy[6]**2))
plt.axis('equal')
plt.show()
