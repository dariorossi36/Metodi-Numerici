# PRECEDENTI, UTILE A RICAVARE SPLITTING ENERGIA TRA PRIMO ECCITATO E FOND E TRA SECONDO ECCITATO E FONDAMENTALE
import pylab
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import string
import matplotlib as mpl

# Cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/Splitting/Bootstrap/"


# Nomi dei file
file_ck = root_dir+'bootstrap_ck.txt'
file_campi = root_dir + 'bootstrap_campi.txt'

# Array in cui salvare i parametri di best fit
E_Y = []
dE_Y = []
E_Y2 = []
dE_Y2 = []

# Valore di N*eta
Neta = 30

k, eta_ck, N_ck, ck, dck, ck2, dck2 = np.loadtxt(file_ck, unpack=True)
eta, N, Y, dY, Y2, dY2 = np.loadtxt(file_campi, unpack=True)

# Prendo solo i primi k_max valori di k e il massimo
k = k[:20]
max_k = int(np.max(k)+1)



# eta_ck e N_ck sono troppo lunghi, perchè ripetono kmax volte lo stesso valore, dato che sono ordinati uguali, uso quindi i valori di eta ed N letti dal campo

# Dichiarazione delle figure

fig1, (ax1, ax3) = plt.subplots(1, 2)
fig2, (ax2, ax4) = plt.subplots(1, 2)

#fig1, ax1 = plt.subplots()
#fig2, ax2 = plt.subplots()

# Separo i c(k) in base all'eta a cui fanno riferimento
for j in range(0, len(eta)):
    
    print(rf"Analizzo $\eta$={eta[j]}")
    
    Ckappa = ck[j*max_k:(j+1)*max_k]
    DCkappa = dck[j*max_k:(j+1)*max_k]
    Ckappa2 = ck2[j*max_k:(j+1)*max_k]
    DCkappa2 = dck2[j*max_k:(j+1)*max_k]

    kappa = k*eta[j]
   
   
   
    ## CALCOLO DELLA FUNZIONE DI CORRELAZIONE CONNESSA
    Ckappa2 = Ckappa2 - (Y2[j])**2
    DCkappa2 = np.sqrt((DCkappa2)**2 + (dY2[j])**2)
    
    # Dummy array per disegnare la funzione di plot
    xx = np.linspace(min(kappa), max(kappa), 1000)
    
    # Funzione di fit
    def funz(x,a,b):
        return a*np.exp(-x*b)
        
    ## FIT DEI CK
    
    print("Analizzo C(k) per il campo")
    
    # Valori iniziali per il Ck del campo
    init = (2,1/1.5)
    
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
    

    # salvataggio dei parametri di best fit
    E_Y = np.append(E_Y, popt[1])
    dE_Y = np.append(dE_Y, err[1])
    
    
    ## FIGURA PER IL CAMPO
    

    # Plot dei C(k) per il campo
    
    #ax1.set_title(fr'Funzione di correlazione a due punti connessa per y' + '\n'+fr' per $\eta N$ = {Neta}')
    ax1.set_xlabel(r'$\omega \tau$')
    ax1.set_ylabel(r'$ \langle y(\tau)y(0) \rangle_C$')
    ax1.errorbar(kappa, Ckappa, DCkappa, linestyle = '', marker ='.', label = fr'$\eta$ = {eta[j]: .3f}')
    ax1.plot(xx, funz(xx, *popt))
    ax1.grid()

    ax1.legend()
    ax1.minorticks_on()
    ax1.grid()
    

    ## FIT DEL CK PER IL CAMPO AL QUADRATO
    print("Analizzo C(k) per il campo al quadrato")
    # Valori iniziali per il Ck del campo
    init = (5,1/5)
    
    popt, pcov=curve_fit(funz, kappa, Ckappa2, init, sigma=DCkappa2, maxfev = 1123456)
    
    
    ndof=len(kappa)-len(init)
    chi2=(((Ckappa2-funz(kappa, *popt))/DCkappa2)**2).sum()
    print('Passo zero')
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))
    dxy=DCkappa2

    dx = np.zeros(len(kappa))
    for i in range(0, 3, 1):
        dxy=np.sqrt(dxy**2+dx**2)
        popt, pcov=curve_fit(funz, kappa, Ckappa2, popt, sigma=DCkappa2)
        chi2=(((Ckappa2-funz(kappa, *popt))/dxy)**2).sum()
        print('Passo %d' % i)
        print('popt:', popt)
        err = np.sqrt(pcov.diagonal())
        print('dpopt:', err)
        print('chi2=%f, \nndof=%f' %(chi2, ndof))
    
    
    print('\n\n\n')
    

    
    # salvataggio dei parametri di best fit
    E_Y2 = np.append(E_Y2, popt[1])
    dE_Y2 = np.append(dE_Y2, err[1])
    

    # Plot dei C(k) per il campo al quadrato

    
    #ax2.set_title(fr'Funzione di correlazione a due punti connessa per $y^2$' + '\n'fr' per $\eta N$ = {Neta}')
    ax2.set_xlabel(r'$\omega \tau$')
    ax2.set_ylabel(r'$ \langle y^2(\tau)y^2(0) \rangle_C$')
    ax2.errorbar(kappa, Ckappa2, DCkappa2, linestyle = '', marker ='.', label = fr'$\eta$ = {eta[j]: .3f}')
    ax2.plot(xx, funz(xx, *popt))
    ax2.grid()

    ax2.legend()
    ax2.minorticks_on()
    ax2.grid()

## FIT LINEARE DEI PARAMETRI DI BEST FIT

# Funzione di fit
def gunz(x,m,q):
    return m+q*x

x = np.square(eta)

# Dummy array per disegnare la funzione di plot
xx = np.linspace(min(x), max(x), 1000)

for i in range(2):

    if (i==0):
        y = E_Y
        dy = dE_Y
        label = r'$\frac{E_1 - E_0}{\hbar \omega}$'
        title = 'Splitting tra fondamentale e primo eccitato'
        ax = ax3
    else:
        y = E_Y2
        dy = dE_Y2
        label = r'$\frac{E_2 - E_0}{\hbar \omega}$'
        title = 'Splitting tra fondamentale e secondo eccitato'
        ax = ax4
    


    # Valori iniziali
    init = (2,-0.1)
    
    # Ciclo per minimizzare il chi quadro
    popt, pcov=curve_fit(gunz, x, y, init, dy)
    
    
    ndof=len(x)-len(init)
    chi2=(((y-gunz(x, *popt))/dy)**2).sum()
    print('Passo zero')
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))
    dxy=dy
    
    dx = np.zeros(len(x))
    for i in range(0, 3, 1):
        dxy=np.sqrt(dxy**2+dx**2)
        popt, pcov=curve_fit(gunz, x, y, popt, dxy)
        chi2=(((y-gunz(x, *popt))/dxy)**2).sum()
        print('Passo %d' % i)
        print('popt:', popt)
        err = np.sqrt(pcov.diagonal())
        print('dpopt:', err)
        print('chi2=%f, \nndof=%f' %(chi2, ndof))
    
    
    print('\n\n\n')
    
    
    # Figura

   # ax.set_title(title)
    ax.set_xlabel(r'$\eta^2$')
    ax.set_ylabel(label)
    
    ax.plot(xx, gunz(xx, *popt), color='red')
    ax.errorbar(x, y, dy, marker ='.', linestyle = '')
    ax.minorticks_on()
    ax.grid()
   
plt.show()