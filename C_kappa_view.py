# Cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"
# A program to plot the C(k) for different beta values for a given lattice lenght

import pylab
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from cycler import cycler
from scipy.optimize import curve_fit


#Scelta della dimensione del reticolo da analizzare
#Nlatt = [10,20,30,40,50,60,70]
Nlatt = [30,50,70]

# Scelta della quantità da analizzare (energia o magnetizzazione)
stringa = 'ene'

# Scelta del beta critico
beta_crit = 0.44

# Set the default color cycle
#custom_cycler = (cycler(color=["#E50000","#E50000","#F97306","#F97306", '#9ACD32', '#9ACD32',"#008000", "#008000","#069AF3","#069AF3", "#12239E","#12239E", "#B66DFF","#B66DFF", "#FF99FF","#FF99FF"]))
custom_cycler = (cycler(color=["#E50000","#F97306", '#9ACD32', "#008000","#069AF3","#12239E", "#B66DFF", "#FF99FF"]))
################################################################################################
######################################## FIT DEI CK ############################################
################################################################################################

# Array in cui salvare i parametri di best fit
tau = []
dtau = []

# Figura complessiva per C al variare di Nlatt, al beta critico
fig1, ax1 = plt.subplots(1,1)
ax1.set_title(fr'$\beta={beta_crit}$', size=13, fontweight='bold')
ax1.set_xlabel(r'$\tau$')
ax1.set_ylabel(r'$C(\tau)$')
ax1.set_prop_cycle(custom_cycler)

for n in range(len(Nlatt)):
    # assign directory dei dati
    directory = root_dir + f"Nlatt={Nlatt[n]}/CK_noidec"

    # Figura unica per l'Nlatt specifico
    fig, ax = plt.subplots(1,1)
    ax.set_title(f'Nlatt={Nlatt[n]}', size=13, fontweight='bold')
    ax.set_xlabel('passi di updating')
    ax.set_ylabel(r'$C(k)$')
    ax.set_prop_cycle(custom_cycler)

    # iterate over files in that directory
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            if f.find('txt')!=-1:

                # salvataggio valore di beta
                inizio = f.find(stringa+'(beta=')    # per importare solo i file nuovi
                if(inizio!=-1):

                # nella prima colonna ho il numero di k e nella seconda il c(k) corrispondente e nella terza l'errore
                    Ckappa, DCkappa, kappa =np.loadtxt(f, unpack=True)

                    inizio = inizio + 9
                    fine = inizio + 5

                    beta=float(f[inizio:fine])

                    # fit
                    # Funzione di fit
                    def funz(x,a,b,c):
                        return a*np.exp(-x/b)+c
                    # Valori iniziali per il Ck
                    init = (1,1/3,0)

                    # Ciclo per minimizzare il chi quadro

                    popt, pcov=curve_fit(funz, kappa, Ckappa, init, sigma=DCkappa)


                    ndof=len(kappa)-len(init)
                    chi2=(((Ckappa-funz(kappa, *popt))/DCkappa)**2).sum()
                    print('Passo zero')
                    print('popt:', popt)
                    print('dpopt:', np.sqrt(pcov.diagonal()))
                    print('chi2=%f, \nndof=%f' %(chi2, ndof))
                    dxy=DCkappa

                    dx = np.zeros(len(kappa))
                    for i in range(0, 3, 1):
                        dxy=np.sqrt(dxy**2+dx**2)
                        popt, pcov=curve_fit(funz, kappa, Ckappa, popt, sigma=DCkappa)
                        chi2=(((Ckappa-funz(kappa, *popt))/dxy)**2).sum()
                        print('Passo %d' % i)
                        print('popt:', popt)
                        err = np.sqrt(pcov.diagonal())
                        print('dpopt:', err)
                        print('chi2=%f, \nndof=%f' %(chi2, ndof))


                    print('\n\n\n')


                    # Impostazioni sui tick
                    ax.minorticks_on()

                    ax.grid(color='lightgray')

                    # grafico

                    #dummy array per plot del fit
                    xx = np.linspace(min(kappa), max(kappa), 1000)

                    ax.errorbar(kappa, Ckappa, DCkappa, linestyle = '', marker ='.', markersize = '3', label = fr'$\beta$ = {beta}')
                    #ax.plot(xx, funz(xx, *popt))
                    print(f'beta = {beta}')
                    ax.legend(loc = 'best')

                    # Impostazioni sui tick
                    ax.minorticks_on()

                    ax.grid(color='lightgray')       #linestyle='--', linewidth=0.1)


                    # salvataggio dei parametri di best fit
                    if(beta==beta_crit):
                        tau = np.append(tau, popt[1])
                        dtau = np.append(dtau, err[1])

                        ax1.errorbar(kappa, Ckappa, DCkappa, linestyle = '', marker ='.', markersize = '2', label = r'$N_{latt}$'+fr'= {Nlatt[n]}')
                        ax1.plot(xx, funz(xx, *popt))
                        # prova=(1, 125 , 0)
                        # ax1.plot(xx,funz(xx,*prova))
                        ax1.legend(loc = 'best')


################################################################################################
####################################### ANALISI DEL TAU ########################################
################################################################################################
# Facciamo analisi del tau al beta critico
beta = 0.441

# Funzione di fit
def g(ics,b,d):
    return ics**b+d
# Valori iniziali per il Ck
init2 = (0.5,10)

# Ciclo per minimizzare il chi quadro

popt, pcov=curve_fit(g, Nlatt, tau, init2, sigma=dtau, maxfev=100)


ndof=len(tau)-len(init2)
chi2=(((tau-g(Nlatt, *popt))/dtau)**2).sum()
print('Passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=dtau

# dx = np.zeros(len(tau))
# for i in range(0, 1, 1):
#     dxy=np.sqrt(dxy**2+dx**2)
#     popt, pcov=curve_fit(g, tau, dtau, popt, sigma=dxy, maxfev=10000000)
#     chi2=(((tau-g(Nlatt, *popt))/dxy)**2).sum()
#     print('Passo %d' % i)
#     print('popt:', popt)
#     err = np.sqrt(pcov.diagonal())
#     print('dpopt:', err)
#     print('chi2=%f, \nndof=%f' %(chi2, ndof))


#dummy array per plot del fit
xx = np.linspace(min(Nlatt), max(Nlatt), 1000)
# Grafico dei tau in funzione degli Nlatt
plt.figure()
#plt.title(r'$\tau_{lat}$ al punto critico in funzione di $N_{lat}$')
plt.xlabel(r'$N_{latt}$')
plt.ylabel(r'$\tau$')
plt.xscale("log")
plt.yscale("log")

# Plot dei C(k)
plt.errorbar(Nlatt, tau, marker ='.', linestyle = '')
plt.minorticks_on()

# Plot della funzione di fit
plt.plot(xx, g(xx, *popt))

plt.show()
