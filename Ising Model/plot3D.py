from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pylab
import os
import matplotlib.gridspec as gridspec

# pacchetto per nuovi colori barra sfumata
import cmasher as cmr
# Access colormap through CMasher or MPL
#cmap = cmr.wildfire                   # CMasher
cmap = plt.get_cmap('cmr.tropical')   # MPL


N=50  #lato reticolo

#beta che mi interessano e quanti sono
beta_inte = [0.2, 0.35, 0.43, 0.44, 0.5]    # questi servono per il caso 2D
quanti_beta = 5 # questi servono per il caso 2D

# beta_inte = [0.2, 0.44, 0.5]    # questi servono per il caso 3D
# quanti_beta = 3 # questi servono per il caso 3D

#cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"
dati_dir = root_dir + f"/Nlatt={N}/Lattice/"

# Dichiarazione figura complessiva nel caso 2D
fig = plt.figure(figsize=(12, 9))
gs = gridspec.GridSpec(2, 6)

# Dichiarazione figura complessiva nel caso 3D
#fig1 = plt.figure(figsize=(12, 4))

# Liste indici degli assi
list_i = [0,1,2,3,4,5,6]    # indici degli assi per i subplot della figura 2d
list_j = [1,2,3,0,0]    # indici degli assi per i subplot della figura 3d

# scorre file nella cartella lattice (un file per ogni beta) e contemporaneamente scorre fli indici da usare sugli assi
for i,j,k in zip(list_i, list_j, np.arange(quanti_beta)):

    #nome file
    f = dati_dir + f"Lattice(beta{beta_inte[k]:.3f}).txt"

    # carico file con valore degli spin
    spin=np.loadtxt(f)
    # DICO QUALE VALORE DI BETA STO ANALIZZANDO
    print(f"\n Analizzo il beta: {beta_inte[k]}")
    # CONTEGGIO DI UP E DOWN
    # Get the unique values and their counts
    unique_values, counts = np.unique(spin.flatten(), return_counts=True)

    # Print the results
    for value, count in zip(unique_values, counts):
        #print(f"{value} occurs {count} times")
        print(f"percentuale di {value}: {((count*100)/(N*N))} %")

# FIGURA 2D

    if i < 3:
        ax = plt.subplot(gs[0, 2 * i:2 * i + 2])
    else:
        ax = plt.subplot(gs[1, 2 * i - 5:2 * i + 2 - 5])

    x=np.arange(N)
    y=np.arange(N)
    X,Y=np.meshgrid(x,y)
    #ax.scatter(X,Y,spin,cmap = 'coolwarm', c = spin, marker = '^')
    im = ax.imshow(spin,cmap = cmap) # disegno 2D con sfumature
    ax.set_title(fr"$\beta= {beta_inte[k]}$"+" \n"+f" spin up: {((count*100)/(N*N))}%")

# FIGURA 3D
   #  # Aggiungo asse 3D
   #  ax1 = fig1.add_subplot(1,3,j, projection='3d')
   # # Punti 3D scatterati
   #  x=np.arange(N)
   #  y=np.arange(N)
   #  X,Y=np.meshgrid(x,y)
   #
   #  ax1.scatter(X,Y,spin,cmap = 'bwr_r', c = spin, marker = '^')
   #
   #  # Assi e titolo
   #  ax1.set_title(fr"$\beta= {beta_inte[k]}$"+" \n"+f" spin up: {((count*100)/(N*N))}%")
   #  ax1.set_zlabel('Valore dello spin')

# FIGURA 2D - barra dei colori
cbar_ax = fig.add_axes([0.9, 0.045, 0.02, 0.4])    # disegno 2D con sfumature
fig.colorbar(im, cax=cbar_ax, label='Valore dello spin')       # disegno 2D con sfumature


plt.tight_layout()
plt.show()