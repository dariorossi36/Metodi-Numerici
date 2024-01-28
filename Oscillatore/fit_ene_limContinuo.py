# PROGRAMMA PER FARE PLOT E FIT DELL'ENERGIA TOTALE IN FUNZIONE DI 1/BETA*OMEGA = 1/n*eta
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/N_ETA_fisso/Bootstrap/"

# Liste da usare
N_eta = []
ene_continuo = []
Dene_continuo = []

for filename in os.listdir(root_dir):
    f = os.path.join(root_dir, filename)

    if f.find('boot')==-1:
        continue
    else:

        y2 , err_y2 , dy2, err_dy2 , y , err_y , ene , err_ene , eta , Nlatt =np.loadtxt(f, unpack=True)

        # Ciclo per salvare i valori di N anche se hanno numero di cifre diverse
        inizio = f.find('NEta=') + 5
        fine = inizio + 5
        NEta = ''
        string = f[inizio:fine]
        for element in string:
            if (element.isdigit() == True):         # controlla se l'elemento nella stringa è un numero
                NEta = NEta + element

        NEta=float(NEta)
        N_eta = np.append(N_eta, NEta)
        print(NEta)
        #Normalizzazione energia
        ene = ene + 1/(2*eta)

        # Funzione di fit
        def f(x,a,b):
            return a+b*x**2

        # Valori iniziali
        init = (1,1)


        # Ciclo per minimizzare il ci quadro
        popt, pcov=curve_fit(f, eta, ene, init, sigma=err_ene)


        ndof=len(eta)-len(init)
        chi2=(((ene-f(eta, *popt))/err_ene)**2).sum()
        print('Passo zero')
        print('popt:', popt)
        print('dpopt:', np.sqrt(pcov.diagonal()))
        print('chi2=%f, \nndof=%f' %(chi2, ndof))
        dxy=err_ene
        dx = np.zeros(len(eta))
        for i in range(0, 3, 1):
            dxy=np.sqrt(dxy**2+dx**2)
            popt, pcov=curve_fit(f, eta, ene, popt, sigma=dxy)
            chi2=(((ene-f(eta, *popt))/dxy)**2).sum()
            print('Passo %d' % i)
            print('popt:', popt)
            print('dpopt:', np.sqrt(pcov.diagonal()))
            print('chi2=%f, \nndof=%f' %(chi2, ndof))


        print('\n\n\n')

        ene_continuo=np.append(ene_continuo,popt[0])
        Dene_continuo=np.append(Dene_continuo,np.sqrt(pcov.diagonal())[0])


       ## Grafico dell'energia potenziale bootstrappata in funzione del valore di eta

        # Dummy array per disegnare la funzione di plot
        xx = np.linspace(min(eta), max(eta), 1000)


        # Figura
        plt.figure()
        plt.title(fr'Energia totale, $N\eta$={NEta}')
        plt.xlabel(r'$\eta$')
        plt.ylabel(r'E totale')

        plt.plot( xx, f(xx, *popt), color='red')
        plt.errorbar(eta, ene, err_ene, marker ='.', linestyle = '')
        plt.minorticks_on()


##  PLOT DELL'ENERGIA IN FUNZIONE DI 1/N eta E FIT
# Funzione di fit
def f(x,a, b):
    return a + b/(np.exp(1/x)-1)

# Valori iniziali
init = (0.5, 1)

ene_continuo= ene_continuo+0.0005
# Ciclo per minimizzare il chi quadro
popt, pcov=curve_fit(f, 1/N_eta, ene_continuo, init, sigma=Dene_continuo )


ndof=len(N_eta)-len(init)
chi2=(((ene_continuo-f(1/N_eta, *popt))/Dene_continuo)**2).sum()
print('Passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=Dene_continuo
dx = np.zeros(len(N_eta))
for i in range(0, 3, 1):
    dxy=np.sqrt(dxy**2+dx**2)
    popt, pcov=curve_fit(f, 1/N_eta, ene_continuo, popt, dxy)
    chi2=(((ene_continuo-f(1/N_eta, *popt))/dxy)**2).sum()
    print('Passo %d' % i)
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))


print('\n\n\n')


## Grafico dell'energia potenziale bootstrappata in funzione del valore di eta

# Dummy array per disegnare la funzione di plot
xx = np.linspace(min(N_eta), max(N_eta), 1000)


# Figura
fig, (ax,ax1)=plt.subplots(2,1, figsize=(18,9), height_ratios=[2,1])

ax.set_ylabel(r'E/$\hbar \omega$',fontsize=13)

ax.plot( 1/xx, f(1/xx, *popt), color='red')
ax.errorbar(1/N_eta, ene_continuo, Dene_continuo, marker ='.', linestyle = '')
ax.minorticks_on()

# Residui normalizzati


r = (ene_continuo-f(1/N_eta,*popt))/dxy
ax1.errorbar( 1/N_eta,r, linestyle='', marker='.')
#ax1.title('Residui normalizzati')
ax1.set_xlabel(r'$\frac{1}{\eta N} $',fontsize=20)
ax1.set_ylabel('Residui norm',fontsize=13)

plt.show()
