import numpy as np
import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

prefix = "/media/user/9EAEE48CAEE45DF1/cyl_temp/p-bed-data/1e4/"
prefix = "../../post/"
prefix = "../"
prefix = prefix + "test_swi1_"
num = np.arange(0,9+1)
Q = np.zeros(len(num))+np.nan
nn = np.zeros(len(num)) + np.nan
for i in range(len(num)):
	nn[i] = i
	name = prefix + str(num[i]).zfill(4) + ".h5"
	print("start process ", name)
	f = h5.File(name)
	Nx = int(f['Nx'][0])
	Ny = int(f['Ny'][0])
	V = np.array(f['Velocity_0'])
	vx = V[0:-2:3]
	vy = V[1:-1:3]
	vx = vx.reshape((Ny,Nx))
	vvx = np.average(vx,axis=1)
	plt.plot(vvx,np.arange(Ny))
	plt.savefig('V.png',dpi=500)
	plt.clf()

	Q[i] = np.sum(vvx)
	plt.plot(nn,Q,'*')
	name = str(i)+'.png'
	plt.savefig(name,dpi=500)
	plt.clf()
