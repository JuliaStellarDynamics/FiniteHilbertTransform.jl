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

# open the file
f = h5py.File('data/data_chebyshev_Plasma_Ku_25_qSELF_0.5_xmax_20.0.hf5','r')
#f = h5py.File('data/data_legendre_Plasma_Ku_25_qSELF_0.5_xmax_20.0.hf5','r')
print(list(f.keys()))

Xi = np.sqrt(np.array(f['tabIminusXi_real'][()])**2.+(1.j*np.array(f['tabIminusXi_imag'][()])**2.)).reshape(801,300)
R = np.array(f['tabomega_real'][()]).reshape(801,300)
I = np.array(f['tabomega_imag'][()]).reshape(801,300)


plt.figure(figsize=(4,4))

plt.contour(R,I,np.log10(1.-Xi),36,colors='k')

plt.xlabel('Re[$\omega$]')
plt.ylabel('Im[$\omega$]')
plt.tight_layout()
plt.savefig('figures/plasma_demo.png')
