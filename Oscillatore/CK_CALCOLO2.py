# PROGRAMMA PER CALCOLARE LA FUNZIONE DI CORRELAZIONE A DUE PUNTI e VALOR MEDIO DEL CAMPO E DEL CAMPO AL QUADRATO
# Restituisce dei file diversi per ogni k e per ogni eta

import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os


# Definizione della cartella principale
numero = 30
scelta = [0.2, 0.4, 0.6, 0.8]
root_dir = f"/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/"
dati_dir = root_dir + "Splitting/Campo/C(k)/"


# Definisco la lunghezza dell'array dei k
max_k = 20


scelta=[]

# Scorro i file dei dati bootstrappati
for filename in sorted(os.listdir(dati_dir)):
    f = os.path.join(dati_dir, filename)

    if f.find('<Y>')!=-1:

        # Salvare valore di eta corrente
        inizio = f.find('eta=') + 5
        fine = inizio + 4
        eta=float(f[inizio:fine])
        scelta=np.append(scelta,eta)
        

print(scelta)

# Scorro i file dei dati bootstrappati
for eta in scelta:
    ipsilon= dati_dir + f"<Y>(eta={eta:.3f},N_Eta={numero}).txt"
    ipsilon2= dati_dir + f"<Y2>(eta={eta:.3f},N_Eta={numero}).txt"
    
    Y=np.loadtxt(ipsilon,unpack=True)
    Y2=np.loadtxt(ipsilon2,unpack=True)
    
    for k in range(max_k):
        f2= dati_dir + f"Ck_Y2(eta={eta:.3f},N_Eta={numero},k={k}).txt"
        f= dati_dir + f"Ck_Y(eta={eta:.3f},N_Eta={numero},k={k}).txt"

        ck=np.loadtxt(f,unpack=True)
        ck2=np.loadtxt(f2,unpack=True)
        
        CK=ck-Y
        CK2=ck2-Y
            

    # Salvataggio dati in un file txt
        nomefile = dati_dir + f"Ck_Y_nuovo(eta={eta:.3f},N_Eta={numero},k={k}).txt"

        nomefile2 = dati_dir + f"Ck_Y2_nuovo(eta={eta:.3f},N_Eta={numero},k={k}).txt"

        np.savetxt(nomefile, CK, delimiter="\t")
        np.savetxt(nomefile, CK2, delimiter="\t")