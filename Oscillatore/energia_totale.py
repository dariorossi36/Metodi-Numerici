# PROGRAMMA PER FARE PLOT E FIT DELL'ENERGIA TOTALE IN FUNZIONE DI 1/BETA*OMEGA = 1/n*eta
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
from uncertainties import ufloat

# Definizione della cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/N_variabile/"
cartel = "Bootstrap/"

directory = root_dir + cartel
# # Definizione degli array da utlizzare
# Nlatt = []
# ene = []
# err_ene = []

# ## Inizializzazione array
# y2 = []
# err_y2 = []
# Nlatt = []
# y2 = []
# dy2 = []
# err_y2 = []
# err_dy2 = []


# # ## iterate over files in that directory
# # for filename in os.listdir(directory):
# #     f = os.path.join(directory, filename)
# #     if f.find('X2')==-1:
# #         continue
# #     else:
# #         # i resampling sono scritti per riga. Leggo il valore di tutti i resampling del campo (è una matrice)
# #         x=np.loadtxt(f, unpack=False)
# #         print(f)

# #         # salvataggio del valore di N a partire dalla stringa
# #         inizio = f.find('N=') + 2
# #         fine = inizio + 4
        
# #         # Ciclo per salvare i valori di N anche se hanno numero di cifre diverse
# #         N = ''
# #         string = f[inizio:fine]
# #         for element in string:
# #             if (element.isdigit() == True):         # controlla se l'elemento nella stringa è un numero
# #                 N = N + element
   
# #         N=float(N)
# #         Nlatt = np.append(Nlatt, N)

# #         # salvataggio del valore di N a partire dalla stringa
# #         inizio = f.find('eta=') + 4
# #         fine = inizio +5

# #         eta=float(f[inizio:fine])

# #         # Analisi della variabile bootstrappata
# #         m_y2=np.mean(x,axis=1) # array con la media delle energie per ogni ricampionamento a eta fisso
# #         y2_media=np.mean(m_y2)     # media delle medie dei ricampionamenti a eta fisso (ovvero media di tutta la matrice)
# #         y2=np.append(y2,y2_media)     # array delle medie cambiando eta

# #         erre = np.std(m_y2, ddof=1)    # errore sulla media di tutti i ricampionamenti a eta fisso
# #         err_y2 = np.append(err_y2, erre)      # errore sulle medie al variare degli eta

# # ## Conversione dlele liste ottenute in array di numpy

# # dy2 = np.array(y2)
# # err_dy2 = np.array(err_y2)
# # Nlatt = np.array(Nlatt)

# # ## iterate over files in that directory
# # for filename in os.listdir(directory):
# #     f = os.path.join(directory, filename)
# #     if f.find('DX2')==-1:
# #         continue
# #     else:
# #         # i resampling sono scritti per riga. Leggo il valore di tutti i resampling del campo (è una matrice)
# #         x=np.loadtxt(f, unpack=False)
# #         print(f)

# #         # Analisi della variabile bootstrappata
# #         m_dy2=np.mean(x,axis=1) # array con la media delle energie per ogni ricampionamento a eta fisso
# #         dy2_media=np.mean(m_dy2)     # media delle medie dei ricampionamenti a eta fisso (ovvero media di tutta la matrice)
# #         dy2=np.append(dy2,dy2_media)     # array delle medie cambiando eta

# #         erre = np.std(m_dy2, ddof=1)    # errore sulla media di tutti i ricampionamenti a eta fisso
# #         print(erre)

# #         err_dy2 = np.append(err_dy2, erre)      # errore sulle medie al variare degli eta

# # ## Conversione dlele liste ottenute in array di numpy
# # dy2 = np.array(y2)
# # err_dy2 = np.array(err_y2)



# # # # Plot del termine cinetico rinormalizzato
# # # partenza=2
# # # y2=y2[partenza:]
# # # eta=eta[partenza:]
# # # #print(f'Nlatt {Nlatt}')
# # # Nlatt=Nlatt[partenza:]
# # # err_y2=err_y2[partenza:]

# # dy2_norm = -dy2/(2*eta**2) + 1/(2*eta)
# # y2=y2/2
# # ene_tot=np.add(dy2_norm,y2)
# # print(len(ene_tot))
# # #print(Nlatt)
# # ascissa=(eta*Nlatt)
# # print(len(ascissa))

# # # Figura
# # plt.figure(1)
# # plt.title('OSCILLATORE ARMONICO \n Termine cinetico')
# # plt.xlabel(r'$\eta$')
# # plt.ylabel(r'$\frac{<\Delta y^2>}{2}$')

# # #plt.plot( xx, f(xx, *popt), color='red')
# # plt.errorbar(ascissa, ene_tot, marker ='.', linestyle = '')
# # plt.minorticks_on()

# # Ciclo per scorrere i file a N diversi (eta è fisso)
# for filename in os.listdir(directory):
#     f = os.path.join(directory, filename)
    
#     if f.find('energy')==-1:
#         continue
#     else:
#         # i resampling sono scritti per riga
#         x=np.loadtxt(f, unpack=False)
#     #faccio deviazione standard delle medie

#     # x=np.loadtxt(r"\\wsl.localhost\Ubuntu\home\dario\Documents\Metodi\Modulo1\Ising\Bootstrap\Nlatt=30\Bootstrappo\Magnetizzazione\resample_binned0.350.txt", unpack=False)
        
#         a=np.mean(x,axis=1)
#         b=np.mean(a)    # dato che plotto dell'energia totale
        
#         ene = np.append(ene, b)

#         err = np.std(a, ddof = 1)
#         err_ene = np.append(err_ene, err)

#         # Ciclo per salvare i valori di N anche se hanno numero di cifre diverse
#         inizio = f.find('N=') + 2
#         fine = inizio + 4
#         N = ''
#         string = f[inizio:fine]
#         for element in string:
#             if (element.isdigit() == True):         # controlla se l'elemento nella stringa è un numero
#                 N = N + element

#         N=float(N)
#         Nlatt = np.append(Nlatt, N)
            
#         inizio = f.find('eta=') + 4
#         fine = inizio + 5
        
#         # Ciclo per salvare i valori di N anche se hanno numero di cifre diverse
#         eta = ''
#         string = f[inizio:fine]
#         for element in string:
#             if (element.isdigit() == True or ('.' in element) == True):         # controlla se l'elemento nella stringa è un numero
#                 eta = eta + element
#                 print(element)
#         eta=float(eta)
        
# ene = np.array(ene)
# err_ene = np.array(err_ene)
# Nlatt_eta = (np.array(Nlatt)*eta)


# #stampo su file tutti i valori che mi interessa plottare
# zippo = zip(Nlatt_eta, ene, err_ene)
# nfile = root_dir + f"/ene_tot.txt"
# zippo=list(zippo)
            
# np.savetxt(nfile, zippo, delimiter="\t", header='N*eta/energia/err_energia')

#############################################################################################################
##############################      ANALISI DATI           ##################################################
#############################################################################################################
nfile = directory + f"/bootstrap_eta=0.008.txt"
y2 , err_y2 , dy2, err_dy2 , y , err_y , ene , err_ene , eta , Nlatt  = np.loadtxt(nfile, unpack=True)


Nlatt_eta=Nlatt*eta
#ene=ene+0.5*Nlatt_eta

## outliers
# np.delete(ene,[3,4,5])
# np.delete(err_ene,[3,4,5])
# np.delete(Nlatt_eta,[3,4,5])

# Funzione di fit
def f(x,a,b):
    return a+1/2+b/(np.exp(x)-1)

# Valori iniziali
init = (-50,1)

# Ciclo per minimizzare il ci quadro
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


## Grafico dell'energia potenziale bootstrappata in funzione del valore di eta

# Dummy array per disegnare la funzione di plot
xx = np.linspace(min(Nlatt_eta), max(Nlatt_eta ), 1000)


# Figura
plt.subplot(211)
plt.title('OSCILLATORE ARMONICO \n Energia totale nel limite al continuo non rinormalizzata')
plt.xlabel(r'$\eta N$')
plt.ylabel(r'U/$\hbar \omega$')

plt.plot( xx, f(xx, *popt), color='red')
plt.errorbar(Nlatt_eta, ene, err_ene, marker ='.', linestyle = '')
plt.minorticks_on()

# Residui normalizzati
plt.subplot(212)
r = (ene-f(Nlatt_eta,*popt))/dxy
plt.errorbar( Nlatt_eta,r, linestyle='', marker='.')
plt.title('Residui normalizzati')
plt.xlabel(r'$\eta N$')
plt.ylabel('(dati - modello)/errore')

plt.show()
