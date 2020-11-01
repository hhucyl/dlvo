import numpy as np 
import h5py as h5
import matplotlib 
matplotlib.use('Agg')
from matplotlib.patches import Ellipse, Circle
import matplotlib.pyplot as plt 
import matplotlib.cm as cm

prefix = 'test_pip_'
d = 2

num = np.arange(0,999+1)

for i in range(len(num)):
    name = prefix + str(num[i]).zfill(4) + ".h5"
    print("start process ", name)
    f = h5.File(name)
    PR = np.array(f['PR'])
    PX = np.array(f['Pposition'])
    px = PX[0:-2:3]
    py = PX[1:-1:3]
    Nx = int(f['Nx'][0])
    Ny = int(f['Ny'][0])
    con = np.array(f['Con'])
    con = con.reshape((Ny,Nx))

    ax = plt.figure().add_subplot(111)
    plt.pcolormesh(con, shading='gouraud',cmap='jet')
    plt.clim(vmin=0, vmax=20)
    for j in range(np.size(PX)/6):
        ppx = px[j]
        ppy = py[j]
        cir = Circle(xy=(ppx,ppy), radius=PR[j]-d, alpha=1.0, facecolor='k')
        ax.add_patch(cir)
    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()
    name = 'c' + str(i) + '.jpg'
    plt.savefig(name, dpi=500)
    plt.clf()
    plt.close()