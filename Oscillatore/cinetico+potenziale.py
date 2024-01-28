import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


## Inizializzazione array
eta = []   
y2 = []   #differenza posizione, generata dalla simulazione
y2_norm = []     #termine cinetico vero e proprio con normalizzazione: 1/2*eta - (deltaY)^2/2*eta^2
err_y2 = []
err_y2_norm = []
Nlatt = []



## assegno cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/"

# Scelta di quale Neta analizzare

NEta = 3


nome_file = root_dir + f'N_ETA_fisso/Bootstrap/bootstrap_NEta={NEta}.txt'

# apertura del file
y2, err_y2, dy2, err_dy2, y, err_y, ene, err_ene, eta, Nlatt  = np.loadtxt(nome_file, unpack = True)


## Grafico del termine cinetico non normalizzato
# Figura
plt.figure(1, figsize=(13,8))
plt.subplot(121)
plt.title(fr'Termine cinetico per $N\eta = {NEta}$')
plt.xlabel(r'$\eta$')
plt.ylabel(r'- $\frac{\langle \Delta y^2 \rangle}{2}$')

plt.errorbar(eta, -dy2/(2*eta**2), err_dy2/(2*eta**2), marker ='.', linestyle = '')
plt.grid()
plt.minorticks_on()

## Funzione di fit (è la stessa per termine cinetico e potenziale)

# Funzione di fit
def f(x,a,b):
    return a*x**2+b

# Valori iniziali
init = (0.65,0.5)


## Fit del termine cinetico
# Definizione del termine cinetico normalizzato
dy2 = - dy2/(2*eta**2) + 1/(2*eta)
err_dy2 = err_dy2/(2*eta**2)

# Funzione di fit
def f(x,a,b):
    return a*x**2+b

# Valori iniziali
init = (0.65,0.5)

# Ciclo per minimizzare il chi quadro
popt, pcov=curve_fit(f, eta, dy2, init, err_dy2)

ndof=len(eta)-len(init)
chi2=(((dy2-f(eta, *popt))/err_dy2)**2).sum()
print('Passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=err_dy2
dx = np.zeros(len(eta))

for i in range(0, 3, 1):
    dxy=np.sqrt(dxy**2+dx**2)
    popt, pcov=curve_fit(f, eta, dy2, popt, dxy)
    chi2=(((dy2-f(eta, *popt))/dxy)**2).sum()
    print('Passo %d' % i)
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))


print('\n\n\n')


# Dummy array per disegnare la funzione di plot
xx = np.linspace(min(eta), max(eta), 1000)

plt.subplot(122)
plt.title(fr'Termine cinetico normalizzato per $N\eta = {NEta}$')
plt.xlabel(r'$\eta$')
plt.ylabel(r'- $\frac{\langle \Delta y^2 \rangle}{2 \eta^2} + \frac{1}{2 \eta}$')

plt.plot( xx, f(xx, *popt), color='red')
plt.errorbar(eta, dy2,  err_dy2, marker ='.', linestyle = '')
plt.minorticks_on()
plt.grid()


## Fit del termine potenziale

# Definizione del termine potenziale
y2 = y2/2
err_y2 = err_y2/2


# Ciclo per minimizzare il chi quadro


popt, pcov=curve_fit(f, eta, y2, init, err_y2)

ndof=len(eta)-len(init)
chi2=(((y2-f(eta, *popt))/err_y2)**2).sum()
print('Passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=err_y2
dx = np.zeros(len(eta))
for i in range(0, 3, 1):
    dxy=np.sqrt(dxy**2+dx**2)
    popt, pcov=curve_fit(f, eta, y2, popt, dxy)
    chi2=(((y2-f(eta, *popt))/dxy)**2).sum()
    print('Passo %d' % i)
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))


print('\n\n\n')


## Grafico dell'energia potenziale bootstrappata in funzione del valore di eta

# Dummy array per disegnare la funzione di plot
xx = np.linspace(min(eta), max(eta), 1000)

# Figura
plt.figure()
plt.title(rf'Termine potenziale per $N\eta = {NEta}$')
plt.xlabel(r'$\eta$')
plt.ylabel(r'$\frac{\langle y^2 \rangle}{2}$')

plt.plot( xx, f(xx, *popt), color='red')
plt.errorbar(eta, y2, err_y2, marker ='.', linestyle = '')
plt.minorticks_on()
plt.grid()


## Grafico con energia cinetica e potenziale sovrapposte
plt.figure()
plt.subplot(211)
plt.title(rf'Differenza termine cinetico e potenziale')
plt.ylabel('Contributi all\'energia')

plt.errorbar(eta, dy2,  err_dy2, marker ='.', linestyle = '', label = 'cinetica')
plt.errorbar(eta, y2, err_y2, marker ='.', linestyle = '', label = 'potenziale')
plt.legend()
plt.grid()
plt.minorticks_on()

plt.subplot(212)
plt.xlabel(r'$\eta$')
plt.ylabel('cinetica-potenziale')
plt.errorbar(eta, dy2-y2, np.sqrt((err_y2)**2+(err_dy2)**2), marker ='.', linestyle = '')
plt.plot( xx, xx*0, color='red')
plt.minorticks_on()
plt.grid()


plt.show()