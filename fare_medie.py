import pylab
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit


# cartella principale
root_dir = fr"/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"

Nlatt=[10,30,50, 70]

for n in range(len(Nlatt)):
    beta=[]
    magne=[]
    err_magne=[]
    # assign directory dei dati
    directory = root_dir + f"Nlatt={Nlatt[n]}/Risultati/"

    for filename in sorted(os.listdir(directory)):
        f = os.path.join(directory, filename)
        if (f.find('txt')!=-1):
            mag, ene, it, b=np.loadtxt(f,unpack=True)

            mag_media=np.mean(mag)
            mag_err=np.std(mag)/len(mag)
            beta=np.append(beta,b[1])
            magne=np.append(magne, mag_media)
            err_magne=np.append(err_magne, mag_err)

    filesaving=root_dir + f"Nlatt={Nlatt[n]}/medi_mag.txt"
    with open(filesaving, 'w') as f:
        for a,b,c in zip(magne, err_magne, beta):
            print("%f  %f  %f" % (a, b, c), file = f)

    f.close()