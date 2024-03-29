# PROGRAMMA VECCHIO PER FARE MEDIA E DEVIAZIONE STANDARD DEI CK BOOTSTRAPPATI DAL ck_bootstrap.c
#Calcola anche la funzione di correlazione a due punti connessa, ovvero la funzione di correlazione normale meno il valor medio di Y o di Y^2 a seconda del caso
import pylab
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

root_dir = f"/home/dario/Documents/UNI/Metodi/Modulo2/Oscillatore/"
dati_dir = root_dir + "Splitting/Bootstrap/"

# Array degli eta scelti
eta_scelto = []

j=0

# Valori di eta
for filename in sorted(os.listdir(dati_dir)):
        f = os.path.join(dati_dir, filename)
        if ((f.find('txt')!=-1) and (f.find('<Y>')!=-1)):
            eta = ''
            inizio1 = f.find('(eta=') + 5
            fine1 = inizio1 + 5
            string1 = f[inizio1:fine1]
            for element in string1:
                if(element.isdigit()==True or ('.' in element) == True):
                    eta= eta + element
                    
            eta = float(eta)
            eta_scelto = np.append(eta_scelto,eta)
print(eta_scelto)
y2=[0.5039, 0.4925, 0.4806, 0.4683, 0.4532]
for i in range(len(eta_scelto)):
    # Array delle correlazioni per campo e campo quadrato
    Ck_Y = []
    dCk_Y = []
    Ck_Y2 = []
    dCk_Y2 = []
    # Array del k
    kappa = []
    # Array del campo o campo quadrato mediato sul cammino (NON SUL TEMPO DEL METROPOLIS!)
    Y = []
    dY = []
    Y2 = []
    dY2 = []
    
    for filename in sorted(os.listdir(dati_dir)):
        f = os.path.join(dati_dir, filename)
        if f.find('txt')!=-1:
            #Salvataggio del valore di eta
            eta = ''
            inizio1 = f.find('(eta=') + 5
            fine1 = inizio1 + 5
            string1 = f[inizio1:fine1]
            for element in string1:
                if(element.isdigit()==True or ('.' in element) == True):
                    eta= eta + element
                    
            eta = float(eta)

            
            if (eta==eta_scelto[i]):
                print(f"sto per aprire x, passo {j} \n")
                j+=1
                   
                x=np.loadtxt(f, unpack=False)
                
                print(np.shape(x))
                
                media = np.mean(x, axis = 1) # Media di un resampling (riga)
                
                if(f.find('ck_(')!=-1):
                    
                    Ck_Y = np.append(Ck_Y, np.mean(media))
                    dCk_Y = np.append(dCk_Y, np.std(media)) # Std della media di ogni resampling

                    k=''
                    inizio = f.find(',k=') + 3
                    fine = inizio + 4
                    string = f[inizio:fine]
                    for element in string:
                        if (element.isdigit() == True):         # controlla se l'elemento nella stringa è un numero
                            k = k + element
                    k = float(k)
                    kappa = np.append(kappa, k)
                    
                if(f.find('ck_2')!=-1):
                    
                    Ck_Y2 = np.append(Ck_Y2, np.mean(media))
                    dCk_Y2 = np.append(dCk_Y2, np.std(media))
                    
                if(f.find('<Y>')!=-1):
                    
                    Y = np.mean(media)
                    dY = np.std(media)
                
                if(f.find('<Y2>')!=-1):
                    
                    Y2 = np.mean(media)
                    dY2 = np.std(media)
                    
    np.asarray(Ck_Y) 
    np.asarray(dCk_Y)
    np.asarray(Ck_Y2)
    np.asarray(dCk_Y2)
                   
    ## CALCOLO DELLA FUNZIONE DI CORRELAIZONE CONNESSA
    Ck_Y = Ck_Y 
    #dCk_Y = np.sqrt(np.square(dCk_Y)+np.square(dY))
    Ck_Y2 = Ck_Y2 - (y2[i])**2
    #dCk_Y2 = np.sqrt(np.square(dCk_Y2)+np.square(dY2)) 
    
    print(f"ck: {len(Ck_Y)}, dck: {len(dCk_Y)}, k: {len(kappa)}")

    # Salvataggio dei dati
    nomefile = root_dir + f"Splitting/CK(eta={eta_scelto[i]}).txt"
    full_array = np.stack([kappa, Ck_Y, dCk_Y, Ck_Y2, dCk_Y2], axis=1) # Unire gli array per salvarli
    np.savetxt(nomefile, full_array, delimiter="\t", header="k\mediaCk_Y\stdCk_Y\mediaCk_Y2\stdCk_Y2" , comments='#')
