# A program that shows the gaussian distribution simulated via the metrogauss.c program.
# It also prints out the acceptance.

import pylab
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes




#y2,Dy2,field,iter=np.loadtxt(fr'C:\Users\aministratore\Documents\Università\Magistrale\Metodi Numerici\Modulo-3\Nuova_run\Oscillo\misure.txt', unpack=True)
Neta = 3
Eta=[1.00, 0.5, 0.1]
# Scelta del numero di bin
binning = 100

fig, axs = plt.subplots(1, 3)
fig.suptitle(r'Istogramma delle posizioni', size=15, weight='bold')


#importo il file del campo genrato dalla simulazione per lo splitting dei livelli energetici da cui prenderò una sola riga in quanto mi interessa 
#il campo per solo un cammino 
for i in range(len(Eta)):
    x=np.loadtxt(fr'/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Campo/Valoricampo(Eta={Eta[i]:.3f})(NEta={Neta}).txt', unpack=False)


    field = x.flatten()

    ## Figura

    #agg = fig.add_subplot(111)   # per aggiungere fig in fig
    axs[i].hist(field, binning, color='lightpink', density=True)
    axs[i].set_title(fr"$\eta={Eta[i]:.1f}$ ")
    axs[i].tick_params(which='minor', length=4, color='r')

    # seleziono i dati che voglio plottare zoomati: per fare questo uso la funzione np.histogram, che mi dice posizione e altezza di ogni bin
    y,x = np.histogram(field, binning, density=True)
    #Taglio i dati così ottenuti, selezionando i bin vicini al picco

    idx = np.argwhere(y>=0.5)

    x_field = (np.take(x, idx)).flatten()
    #tmp = x_field + (x_field[1]-x_field[0])
    #x_field = (tmp+x_field)/2
    y_field = (np.take(y, idx)).flatten()


    xx=np.linspace(-3,3,10000)
    '''
    A=1.7
    B=0
    C=0.05
    '''

    A=1/(np.pi**(1/4))
    B=0
    C=1


    ground=A*np.exp(-((xx+B)**2)/(2*C))
    axs[i].errorbar(xx,abs(ground)**2, marker ='', color='purple', linestyle = '-')


        # FIGURA ZOOMATA
    # location for the zoomed portion
    axins3 = inset_axes(axs[i], width="30%", height=1., loc=1)

    # plot the zoomed portion
    axins3.step(x_field, y_field,color = 'lightpink', where= 'post') #, c = 'k')
    # dummy array tagliato
    dist = (x_field[1]-x_field[0])
    xx_cut = np.linspace(np.min(x_field),np.max(x_field),10000)



    # funzione d'onda com dummy array tagliato
    ground_cut=A*np.exp(-((xx_cut+B)**2)/(2*C))

    axins3.fill_between(x_field, y_field, y2=min(y_field), alpha=0.8, color = 'lightpink', step = 'post') #y[idx[0]-1], abs(ground_cut[0])**2

    axins3.errorbar(xx_cut,abs(ground_cut)**2, marker ='', color='purple', linestyle = '-')


    axs[i].tick_params(which='minor', length=4, color='r')

pylab.show()