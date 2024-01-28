import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


## assegno cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/"
n_file = 'N_ETA_fisso/ene_continuo.txt'

N_eta, a, pot_continuo, a_err, pot_cont_err, chi2, ndof =np.loadtxt(n_file, unpack=True)

plt.figure(1)
plt.errorbar(N_eta,pot_continuo, pot_cont_err, marker ='.', linestyle = '')
plt.minorticks_on()
plt.show()