root_dir = fr"/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"
# A program to plot the means of the rebinned bootstrap

import pylab
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib as mpl
from cycler import cycler
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

Nlatt=[10,20,30,40,50,60,70]

plt.figure()

for n in range(len(Nlatt)):
    beta=[]
    magne=[]
    err_magne=[]
    # assign directory dei dati
    filename = root_dir + f"Nlatt={Nlatt[n]}/medi_mag.txt"
    magne, err_magne, beta=np.loadtxt(filename,unpack=True)    
    
    plt.errorbar(beta,abs(magne),err_magne, linestyle='', marker ='.',label=fr'Nlatt = {Nlatt[n]}')
plt.legend()
plt.show()