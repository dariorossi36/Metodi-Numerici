# PROGRAMMA PER FARE PLOT E FIT DELL'ENERGIA TOTALE IN FUNZIONE DI 1/BETA*OMEGA = 1/n*eta
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

label = 'Continuo/'     # Scelta per limite al continuo ('Continuo') o cambio di temperatura ('Temperatura')

# Cartelle
root_dir ="/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/"


# Inizializzazione della figura

fig, axs = plt.subplots(2,2,height_ratios = [2,1])
fig.suptitle('Energy density in funzione della massa', size = 15, weight='bold')

for i in range(0,2):

    if(i==1):
        relax = 'NoOverRelax/'
        title = 'Senza Over Relaxation'
        print('Senza Over Relaxation')
    else:
        relax = ''
        title = 'Con Over Relaxation'
        print('Con Over Relaxation')

    dati_dir = root_dir + label + relax + "Bootstrap/"


    nfile =  dati_dir + f"/bootstrap.txt"
    epsilon, depsilon, mass, dmass, m, Nt, Nx, Ntm, Nxm = np.loadtxt(dati_dir+"bootstrap.txt",unpack=True)
    ## Fit con una parabola

    #calcolo densità di energia su temperatura al quadrato
    epsilon = epsilon*Nt**2
    depsilon = depsilon*Nt**2


    # outlier
    y=epsilon
    x=m
    dy=depsilon
    taglio=2  #numero di dati che voglio tagliare per m grosso

    ordine=np.argsort(x)
    for j in range (taglio):
        outlier = np.argmax(y)
        x = np.delete(x,outlier)
        y = np.delete(y,outlier)
        dy = np.delete(dy,outlier)


    ordine2=np.argsort(x)
    altri_outlier=ordine2[0],ordine2[3], ordine2[6] #anche i dati per m piccola non mi piacciono
    '''
    x = np.delete(x,altri_outlier)
    y = np.delete(y,altri_outlier)
    dy = np.delete(dy,altri_outlier)
    '''

    # Funzione di fit
    def funz(ics,a,b):
        return a*ics**2+b

    # Valori iniziali
    init = (1,np.pi/6)

    # Ciclo per minimizzare il chi quadro


    popt, pcov=curve_fit(funz, x, y, init, dy)


    ndof=len(x)-len(init)
    chi2=(((y-funz(x, *popt))/dy)**2).sum()
    print('Passo zero')
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))
    dxy=dy

    dx = np.zeros(len(x))
    for j in range(0, 3, 1):
        dxy=np.sqrt(dxy**2+dx**2)
        popt, pcov=curve_fit(funz, x, y, popt, dxy)
        chi2=(((y-funz(x, *popt))/dxy)**2).sum()
        print('Passo %d' % i)
        print('popt:', popt)
        err = np.sqrt(pcov.diagonal())
        print('dpopt:', err)
        print('chi2=%f, \nndof=%f' %(chi2, ndof))


    print('\n\n\n')

    # Plot dell'energy density in funzione della massa

    # Titolo ed etichette per gli assi
    labelx = r'$\hat{m}$'
    labely = r'$\frac{\epsilon}{T^2}$'

    #    Dummy array per plot con la funzione di fit
    xx = np.linspace(0, max(x), 1000)

    axs[0][i].set_title(title, size=8)
    axs[0][0].set_ylabel(labely)


    axs[0][i].errorbar(x, y, dy, marker ='.', linestyle = '', color='green',label='dati')
    massa=np.sort(m)
    epsilon_2=np.take(epsilon,ordine)
    depsilon_2=np.take(depsilon,ordine)
    #axs[0][i].errorbar(massa[len(m)-taglio:], epsilon_2[len(m)-taglio:], depsilon_2[len(m)-taglio:], marker ='.', linestyle = '', label='outliers',color='deepskyblue',)
    axs[0][i].plot(xx, funz(xx, *popt),color = 'orange',label='best-fit')
    axs[0][i].plot(xx, np.pi/6*np.ones(len(xx)), color='red', linestyle='--', label=r'$\pi/6$')
    axs[0][i].minorticks_on()

    axs[0][i].legend()
    axs[0][i].grid(linewidth=0.5)

    # residui normalizzati
    axs[1][i].errorbar(x, (y-funz(x, *popt))/dy, marker = '.', linestyle = '')
    axs[1][i].plot(x,np.zeros(len(x)), linestyle = '--')
    axs[1][i].grid(linewidth=0.5)
    axs[1][0].set_ylabel('Residui norm')
    axs[1][i].set_xlabel(labelx)
    #axs[1][i].minorticks_on()


plt.show()
