root_dir = "/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"

import pylab
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from cycler import cycler

Nlatt=50

dati_dir = root_dir + f"Nlatt={Nlatt}/Bootstrap/"

#liste per il plot
bins=[]
err=[]


# Set the default color cycle
#custom_cycler = (cycler(color=["#000000",  "#E50000","#F97306", '#9ACD32',"#008000","#069AF3", "#12239E","#B66DFF","#FF99FF","#FF9999"]))
#custom_cycler = (cycler(color=["b",  "#E50000","y"]))
# Figura unica
fig, ax = plt.subplots(1,1, figsize=(10,8))
#pylab.title(f'Autocorrelazione della magnetizzazione verso la transizione \n reticolo {Nlatt}x{Nlatt}')
#ax.set_xlabel('binning')
ax.set_xlabel('resampling')
ax.set_ylabel('incertezza')
#ax.set_prop_cycle(custom_cycler)

# iterate over files in that directory
for filename in os.listdir(dati_dir):
    f = os.path.join(dati_dir, filename)
    # checking if it is a file
    if os.path.isfile(f):
        #if f.find('res')!=-1:
        if f.find('bin')!=-1:

            # salvataggio valore del binning
            inizio = f.find('bin=')  
            #inizio = f.find('res=')  
            if(inizio!=-1):

                # Importo dati del bootstrap
                ene, err_ene, mag, err_mag, calore, err_calore, susce, err_susce, binder, err_binder, beta = np.loadtxt(f, unpack=True)
                
                scelta=err_susce  #scelta errore da plottare nel caso di variazione binning
                #scelta=err_mag #scelta errore da plottare nel caso di variazione resampling
                taglio=10
                
                inizio = inizio + 4
                fine = inizio + 5
                
                binning=''
                b=(f[inizio:fine])
                for str in b:
                    if (str.isdigit()==True):
                        binning=binning+str
                binning=int(binning)
                mult=12
                for i in range(1):
                    index=np.argwhere(beta==0.44)
                    ax.errorbar(int(binning),scelta[index+i*mult],marker='o', linestyle='', color="b")
                    if (binning==5):
                        ax.errorbar(int(binning),scelta[index+i*mult],marker='o',label=fr'$\beta$={float(beta[index+mult*i])}', linestyle='', color="b")
ax.legend(loc='upper right')
ax.minorticks_on()
ax.grid()




plt.show()