# Definizione della cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"

#Un programma per vedere la storia montecarlo della magnetizzazione e dell'energia
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#specifiche del reticolo e del beta che voglio plottare
L=np.array([10,30, 50,70])
Nlatt=50
beta_pocopiudic=0.448
beta_c=0.44
beta_grosso=0.53
beta_piccolo=0.32

#taglio misure che voglio plottare
taglio_sx=21000
taglio_dx=23000

#importo dati simulazione e inserisco specifiche della stessa
scelta=[beta_piccolo,beta_c, beta_pocopiudic,beta_grosso]


idec=100



fig, ax=plt.subplots(4,2)

fig.suptitle(r'Storia Monte-Carlo e distribuzione Magnetizzazione', size=20, weight='bold')
for j in range(len(scelta)):
    dati_dir = root_dir + f'Nlatt={Nlatt}/Risultati/misure({scelta[j]:.3f}).txt'

    mag,ene,passo,b=np.loadtxt(dati_dir, unpack=True) #la matrice in questione
    #ax[j][0].set_title(fr'$\beta$={scelta[j]}',size=11)
    #ax[j][0].set_xlabel('Passi')
    ax[j][0].set_ylabel('M')
    ax[j][0].set_ylim(-1,1)
    if (scelta[j]==beta_c):
        ax[j][0].errorbar(passo,mag, linestyle = '-', color='m', marker ='.', markersize=1, label=fr'$\beta$={scelta[j]}')
        x_field = passo[taglio_sx:taglio_dx].flatten()
        y_field = mag[taglio_sx:taglio_dx].flatten()
        ax[j][0].legend(loc='lower right')
        # location for the zoomed portion
        axins3 = inset_axes(ax[j][0], width="30%", height=1., loc=1)

        # plot the zoomed portion
        axins3.step(x_field, y_field,color = 'm', where= 'post') #, c = 'k')
        axins3.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, left=False, labelleft=False, right=True, labelright=True)

    else:
        ax[j][0].errorbar(passo,mag, linestyle = '-', color='m', marker ='.', markersize=1, label=fr'$\beta$={scelta[j]}')
        ax[j][0].legend()
    ax[j][0].minorticks_on()
    ax[j][0].grid()
    

    #ax[j][1].set_title(, size=11)
    ax[len(scelta)-1][1].set_xlabel('M')
    ax[j][1].set_ylabel(r'$\mathcal{P}[M]$')
    ax[j][1].hist(mag,70, color='lightpink', label=fr'$\beta$={scelta[j]}')
    ax[j][1].set_xlim(-1,1)
    #y,x = np.histogram(mag,70)
    #ax[j][1].step(x[:-1],y, where='post', color='m')
    ax[j][1].minorticks_on()
    ax[j][1].grid()
    ax[j][1].legend()

fig2, ax2=plt.subplots(1,len(L),figsize=(10,5))
fig2.suptitle(rf'Distribuzione Magnetizzazione a $\beta$={beta_piccolo} ', size=15, weight='bold')
for i in range(len(L)):
    dati_dir = root_dir + f'Nlatt={L[i]}/Risultati/misure({beta_piccolo:.3f}).txt'

    mag,ene,misure,b=np.loadtxt(dati_dir, unpack=True) #la matrice in questione
    ax2[i].set_title(rf'Nlatt={L[i]}', size=11)
    ax2[i].set_xlabel('M')
    ax2[0].set_ylabel(r'$\mathcal{P}[M]$')
    ax2[i].set_xlim(-1,1)
    ax2[i].set_ylim(0,10)
    ax2[i].grid()
    if (L[i]==10):
        ax2[i].hist(mag, 20, color='purple', density=True)
    else:
        ax2[i].hist(mag, 30, color='purple', density=True)
    ax2[i].minorticks_on()
plt.tight_layout()
plt.show()