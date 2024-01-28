# PROGRAMMA PER FARE PLOT E FIT DELL'ENERGIA TOTALE IN FUNZIONE DI 1/BETA*OMEGA = 1/n*eta
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

label = 'Temperatura/Simmetrico'

# Cartelle
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo3/Scalar2D/"
dati_dir = root_dir + label + "/Bootstrap/"

nfile =  dati_dir + f"/bootstrap_mix.txt"
epsilon, depsilon, mass, dmass, m, Nt, Nx, Ntm, Nxm = np.loadtxt(nfile,unpack=True)
# Sottraggo il valore di epsilon corrispondente a Nt più alto (T più bassa), perché non sappiamo il valore zero dell'energia. Quindi l'ultimo

idx = np.argmax(Nt)
print(idx, Nt[idx])

# Calcolo dell'energia con la somma delle osservabili singole

ene = 0.5*(epsilon - epsilon[idx])*Nt**2
dene = 0.5*np.sqrt(depsilon)*Nt**2

# Cambio variabile per comodità
x1 = Nt
y1 = 0.5*epsilon
dy1 = 0.5*depsilon
#y1 = ene[:stop]
#dy1 = dene[:stop]

## Plot dell'energy density in funzione di Nt

# Titolo ed etichette per gli assi
title = 'Energy density in funzione di Nt'
labelx = r'$Nt$'
labely = r'$\epsilon$'


plt.figure()
#plt.title(title)
plt.xlabel(labelx)
plt.ylabel(labely)

plt.errorbar(Nt*m, y1, dy1, marker ='.', linestyle = '', color = 'black')

plt.minorticks_on()

# Normalizzazione dell'energy density
#epsilon[-2]=epsilon[-2]+0.00003
#epsilon[-4]=epsilon[-4]+0.00003
epsilon = (epsilon - epsilon[np.argmax(Nt)])

# Cambio variabile per comodità
stop = len(epsilon)
x2 = (1/(Nt[:stop]*m))
y2 = 0.5*(epsilon[:stop])*Nt**2
dy2 = 0.5*np.sqrt(depsilon[:stop]**2+depsilon[np.argmax(Nt)]**2)*Nt**2

zippo=zip(x1, y1, dy1, x2, y2, dy2)
zippo=list(zippo)
nfile =  root_dir + label + f"/energy_density.txt"
np.savetxt(nfile, zippo, delimiter='\t',  header='#Nt / density / delta(sx) / (1/Ntm) / density*Nt**2 / Delta(sx)')

######################################################################################
##########################  RUNNA QUESTO SE VUOI I GRAFICI ###########################
#######################################################################################


## Plot dell'energy density in funzione di 1/Nt*m e moltiplicata per Nt**2

# Titolo ed etichette per gli assi
title = r'Energy density in funzione di $\frac{1}{Nt*m}$'
labelx = r'$\frac{1}{Nt*\hat{m}}$'
labely = r'$\frac{\epsilon}{T^2}$'


plt.figure(figsize=(8,6))
#plt.title(title)
plt.xlabel(labelx,fontsize=16)
plt.ylabel(labely,fontsize=16)


plt.errorbar(x2, y2, dy2, marker ='.', linestyle = '')
#plt.errorbar(x, y1, marker ='.', linestyle = '', color = 'green')
xx = np.linspace(min(x2), max(x2), 1000)
plt.plot(xx, np.pi/6*np.ones(len(xx)), linestyle = '-')
plt.minorticks_on()
plt.show()