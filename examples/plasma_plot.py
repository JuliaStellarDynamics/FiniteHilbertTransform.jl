"""
plasma_plot.py
    visualise some of the plasma calculations

"""

import numpy as np

# enable hdf5 reading
import h5py

# plotting basics
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# override some parameters for a 'house style'
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['font.weight'] = 'medium'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 0.75
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.width'] = 0.75
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True


filename = 'data/data_Plasma_Analytical.h5'


f = h5py.File(filename, 'r')

epsilon = np.abs(np.array(f['tabIminusXi_real']).reshape([801,300])+1j*np.array(f['tabIminusXi_imag']).reshape([801,300]))
plt.contour(np.array(f['tabomega_real']).reshape([801,300]),\
             np.array(f['tabomega_imag']).reshape([801,300]),\
             np.log10(epsilon),[-2.,-1.5,-1.,-0.5,-0.25,0.],vmin=-2.,vmax=-1.,colors='k')
plt.colorbar()
plt.xlabel("Re[$\omega$]")
plt.ylabel("Im[$\omega$]")
#plt.title("Isochrone mode test")
plt.savefig('figures/plasma_demo.png')


# open the file
f = h5py.File('data/data_chebyshev_Plasma_Ku_200_qSELF_0.5_xmax_20.0.hf5','r')
#f = h5py.File('data/data_legendre_Plasma_Ku_200_qSELF_0.5_xmax_20.0.hf5','r')
print(list(f.keys()))

epsilon = np.abs(np.array(f['tabIminusXi_real']).reshape([801,300])+1j*np.array(f['tabIminusXi_imag']).reshape([801,300]))

R = np.array(f['tabomega_real'][()]).reshape(801,300)
I = np.array(f['tabomega_imag'][()]).reshape(801,300)


plt.figure(figsize=(4,4))

plt.contour(R,I,\
             np.log10(epsilon),[-2.,-1.5,-1.,-0.5,-0.25,0.],vmin=-2.,vmax=-1.,colors='k')

plt.xlabel('Re[$\omega$]')
plt.ylabel('Im[$\omega$]')
plt.tight_layout()
plt.savefig('figures/plasma_demo.png')
