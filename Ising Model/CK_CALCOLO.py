# Definizione della cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"
# PROGRAMMA PER CALCOLARE LE FUNZIONI DI CORRELAZIONE A DUE PUNTI e il valor medio del campo e del campo al quadrato. RESTITUISCE DEI FILE DIVERSI PER OGNI k E PER OGNI eta
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Dimensione reticolo
Nlatt = 10
print(Nlatt)

# Scelta di beta interessanti
beta_inte = [0.44]


# Cartella dei dati da analizzare
scelta = "Energia/"
dati_dir = root_dir + f"Nlatt={Nlatt}/Bootstrap_noidec/" + scelta

# Definizione della cartella delle misure
mis_dir = root_dir + f"Nlatt={Nlatt}/CK_noidec/"

# Definisco la lunghezza dell'array dei k
inizio =0
max_k = 200

# Scorro i file dei dati (sono dati a beta diversi)
for beta in beta_inte:

    # Apro il file
    filename = dati_dir + f"resample_binned{beta:.3f}.txt"
    x=np.loadtxt(filename, unpack=False)

    righe = np.shape(x)[0]  #Numero dei resampling
    N = np.shape(x)[1]  #Numero dei dati (cioè delle colonne della matrice)


    media = np.mean(x, axis=1)  # media di un resampling + il bootstrap mi serve per avere l'errore

    # Inizializzazione degli array
    stdCk = np.zeros(max_k)
    meanCk = np.zeros(max_k)
    kappa = np.zeros(max_k)

    print(f"media: {len(media)}, shape x: {np.shape(x)}")


# Calcolo del Ck per k=0 per poi normalizzare
    Ck_norm = np.zeros(righe)
    for i in range(righe):
        for j in range(N):
            Ck_norm[i] = (x[i][j] - media[i])*(x[i][j] - media[i]) + Ck_norm[i]

        Ck_norm[i] = Ck_norm[i]/N     #Divido per numero di cose sommmate


    meanCk_norm = np.mean(Ck_norm)     # Sto mediando sui resampling, a k fisso
    stdCk_norm = np.std(Ck_norm)


# Calcolo del Ck
    for k in range(inizio,max_k):
        Ck = np.zeros(righe)

        for i in range(righe):
            for j in range(N-k):
                Ck[i] = (x[i][j] - media[i])*(x[i][j+k] - media[i]) + Ck[i]
                #Ck[i] = (x[i][j])*(x[i][j+k]) + Ck[i]

            Ck[i] = Ck[i]/(N-k)     #Divido per numero di cose sommmate


        kappa[k]=int(k)     # salvo il k in un array da salvare nel file txt
        meanCk[k] = np.mean(Ck)     # Sto mediando sui resampling, a k fisso
        stdCk[k] = np.std(Ck)

    # Normalizzo Ck
    stdCk = np.sqrt((stdCk_norm)**2/(meanCk_norm)**2 + (stdCk[inizio:])**2/(meanCk[inizio:])**2) # questo lo faccio prima perché dopo ridefinisco meanCk
    meanCk = meanCk/meanCk_norm
    stdCk = stdCk*meanCk[inizio:]

    # Salvataggio dati in un file txt
    if(filename.find('Energia')!=-1):
        nomefile = mis_dir + f"CK_ene_surplus(beta={beta:.3f}).txt"
    else:
        nomefile = mis_dir + f"CK_mag(beta={beta:.3f}).txt"
    np.savetxt(nomefile, np.column_stack([meanCk[inizio:], stdCk, kappa[inizio:]]), delimiter="\t")

