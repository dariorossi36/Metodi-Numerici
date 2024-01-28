# PROGRAMMA PER FARE PLOT E FIT DELL'ENERGIA TOTALE IN FUNZIONE DI 1/BETA*OMEGA = 1/n*eta
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Array degli eta che mi interessano
eta_scelto = np.array([0.010, 0.008, 0.005, 0.002])
# Array delle divergenze
div = []
Ddiv = []

# Definizione della cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/N_variabile/"
cartel = "Bootstrap"

directory = root_dir + cartel


for i in range (len(eta_scelto)):

    nfile = directory+ f"/bootstrap_eta={eta_scelto[i]:.3f}.txt"
    y2, err_y2, dy2, err_dy2, y, err_y, ene, err_ene, eta, Nlatt = np.loadtxt(nfile, unpack=True)

    Nlatt_eta=Nlatt*eta


    ## GRAFICO E FIT DELL'ENERGIA IN FUNZIONE DI 1/BETA*OMEGA = 1/n*eta
    # Funzione di fit
    def f(x,a,b):
        return a+1/2+b/(np.exp(x)-1)

    # Valori iniziali
    init = (-50,1)

    # Ciclo per minimizzare il ci quadro
    print(f"eta={eta_scelto[i]}")
    popt, pcov=curve_fit(f, Nlatt_eta, ene, init, err_ene)


    ndof=len(Nlatt_eta)-len(init)
    chi2=(((ene-f(Nlatt_eta, *popt))/err_ene)**2).sum()
    print('Passo zero')
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))
    dxy=err_ene
    dx = np.zeros(len(Nlatt_eta))
    for i in range(0, 3, 1):
        dxy=np.sqrt(dxy**2+dx**2)
        popt, pcov=curve_fit(f, Nlatt_eta, ene, popt, dxy)
        chi2=(((ene-f(Nlatt_eta, *popt))/dxy)**2).sum()
        print('Passo %d' % i)
        print('popt:', popt)
        print('dpopt:', np.sqrt(pcov.diagonal()))
        print('chi2=%f, \nndof=%f' %(chi2, ndof))


    print('\n\n\n')

    # Salvataggio dei valori di best fit in un array
    div = np.append(div, popt[0])
    Ddiv = np.append(Ddiv, np.sqrt(pcov.diagonal())[0])



    ## Grafico dell'energia potenziale bootstrappata in funzione del valore di eta

    # Dummy array per disegnare la funzione di plot
    xx = np.linspace(min(Nlatt_eta), max(Nlatt_eta ), 1000)


    # Figura
    plt.figure()
    plt.subplot(211)
    #plt.title(fr'OSCILLATORE ARMONICO \n Energia totale nel limite al continuo non rinormalizzata $\eta = {eta_scelto[i]}')
    plt.xlabel(r'$\eta$')
    plt.ylabel(r'E/$\hbar \omega$')

    plt.plot( xx, f(xx, *popt), color='red')
    plt.errorbar(Nlatt_eta, ene, err_ene, marker ='.', linestyle = '')
    plt.minorticks_on()

    # Residui normalizzati
    plt.subplot(212)
    r = (ene-f(Nlatt_eta,*popt))/dxy
    plt.errorbar( Nlatt_eta,r, linestyle='', marker='.')
    plt.title('Residui normalizzati')
    plt.xlabel(r'$\frac{1}{\eta N}$')
    plt.ylabel('(dati - modello)/errore')



## FIT E PLOT DI 1/(2eta) IN FUNZIONE DI 1/eta
print("inizia il fit del termine divergente")
def g(x,a,b):
    return a*x+b

# Valori iniziali
init = (0.5,0)

# Ciclo per minimizzare il ci quadro
popt, pcov=curve_fit(g, 1/eta_scelto, div, init, Ddiv)


ndof=len(eta_scelto)-len(init)
chi2=(((div-g(1/eta_scelto, *popt))/Ddiv)**2).sum()
print('Passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=Ddiv
dx = np.zeros(len(eta_scelto))
for i in range(0, 3, 1):
    dxy=np.sqrt(dxy**2+dx**2)
    popt, pcov=curve_fit(g, 1/eta_scelto, div, popt, dxy)
    chi2=(((div-g(1/eta_scelto, *popt))/dxy)**2).sum()
    print('Passo %d' % i)
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))


print('\n\n\n')

## Grafico dell'energia potenziale bootstrappata in funzione del valore di eta

# Dummy array per disegnare la funzione di plot
xx = np.linspace(min(eta_scelto), max(eta_scelto), 1000)


# Figura
plt.figure()
plt.subplot(211)
#plt.title("OSCILLATORE ARMONICO \n Fattore correttivo divergente all'energia cinetica")
plt.xlabel(r'$\frac{1}{\eta}$')
plt.ylabel(r'Correzione')

plt.plot( 1/xx, g(1/xx, *popt), color='red')
plt.errorbar(1/eta_scelto, div, Ddiv, marker ='o', linestyle = '')
plt.minorticks_on()

# Residui normalizzati
plt.subplot(212)
r = (div-g(1/eta_scelto,*popt))/dxy
plt.errorbar(1/eta_scelto, r, linestyle='', marker='o')
plt.title('Res norm')
plt.xlabel(r'$\frac{1}{\eta}$')
plt.ylabel('Res norm')

plt.show()