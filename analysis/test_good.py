import h5py
import sys
import numpy as np

prefix = sys.argv[1]
fistart = int(sys.argv[2])
fiend = int(sys.argv[3])+1
fn = np.arange(fistart, fiend)

for i in fn:
	name = prefix + str(i).zfill(4) + '.h5'
	f = h5py.File(name)	
	Nx = np.array(f['Nx'])
	Ny = np.array(f['Ny']);
	vel = np.array(f['Velocity_0'])
	vx = vel[0:-2:3]
	vx = vx.reshape(Ny[0],Nx[0])
	vvx = np.sum(np.abs(np.average(vx,axis=1)))
	print(name, vvx)
	f.close()
