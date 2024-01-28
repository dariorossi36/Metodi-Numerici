root_dir = fr"/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"
# A program to plot the means of the rebinned bootstrap

import pylab
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib as mpl
from cycler import cycler
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


#Dichiaro gli esponenti critici analiticamente corretti, proprio beta doveva chiamarsi vero?
bbeta = 1/8
gamma = 7/4
beta_C = 0.441

# Rimozione del file dell'esponente critico nu, se esiste già (per non sovrascrivere)
filenu = root_dir + "/Exp_crit/exp_nu.txt"
if os.path.exists(filenu):
    os.remove(filenu)

#Scelta dei colori da usare nel grafico
plt.rc('lines', linewidth=1.2)
plt.rc('axes', prop_cycle=(cycler('color', ["#E50000","#F97306", '#9ACD32', "#008000", "#069AF3", "#12239E", "#B66DFF"])))
custom_cycler = (cycler(color=["#E50000","#E50000","#F97306","#F97306", '#9ACD32', '#9ACD32',"#008000", "#008000","#069AF3","#069AF3", "#12239E","#12239E", "#B66DFF","#B66DFF"]))

#Creazione di una figura unica per mettere a confronto reticoli di dimensioni diverse
fig_c, ax_c = plt.subplots(figsize=(10,7))
fig_m, ax_m = plt.subplots(figsize=(10,7))
fig_s, ax_s = plt.subplots(figsize=(10,7))
fig_b, ax_b = plt.subplots(figsize=(10,7))

#Stesse figure di sopra ma riscalate per il Finite Size Scaling
fig_m_FSS, ax_m_FSS = plt.subplots(figsize=(10,7))
fig_s_FSS, ax_s_FSS = plt.subplots(figsize=(10,7))

# location for the zoomed portion
axins3 = inset_axes(ax_b, width="40%", height=1.5, loc=1)

# Resetto ciclo colori per l'asse del binder (punti e linea dello stesso colore)
ax_b.set_prop_cycle(custom_cycler)
axins3.set_prop_cycle(custom_cycler)

#Inizio ciclo per scorrere gli nlatt

for root, dirs, files in os.walk(root_dir):
    dirs = sorted(dirs)
    for i in range(len(dirs)):

        index = dirs[i].find('Nlatt')
        if index !=-1:
            inizio = index + 6
            fine = index + 8
            Nlatt = (dirs[i])[inizio:fine]
            Nlatt=int(Nlatt)

            if(Nlatt!=80):
                filename = root_dir + dirs[i] + "/Bootstrap/bootstrap.txt"
                # Importo dati del bootstrap
                ene, err_ene, mag, err_mag, calore, err_calore, susce, err_susce, binder, err_binder, beta = np.loadtxt(filename, unpack=True)

                # Grafico del Calore specifico bootstrappato in funzione dei beta
                ax_c.errorbar(beta, calore, err_calore, marker ='.', linestyle = '', label=f'Nlatt = {Nlatt}')
                # Grafico del suscettività bootstrappata in funzione dei beta
                ax_s.errorbar(beta, susce, err_susce, marker ='.', linestyle = '',label=f'Nlatt = {Nlatt}')
                # Grafico della magnetizzazione bootstrappata in funzione dei beta
                ax_m.errorbar(beta, mag, err_mag, linestyle = '', marker ='.',label=f'Nlatt = {Nlatt}')

                ##### FINITE SIZE SCALING ####
                #riscalo suscettività ed errore associato
                ssusc=susce/Nlatt**gamma
                ssusc_err=(err_susce**2)/Nlatt**gamma
                t=(beta-beta_C)*Nlatt

                #riscalo magnetizzazione e errore associato
                mmag=np.abs(mag/Nlatt**(-bbeta))
                err_mmag=err_mag/Nlatt**(-bbeta)
                tm=(beta-beta_C)*Nlatt

                # Grafico del suscettività in Finite State Size in funzione dei beta
                ax_s_FSS.errorbar(t, ssusc, ssusc_err, marker ='.', linestyle = '',label=f'Nlatt = {Nlatt}')
                # Grafico della magnetizzazione in Finite State Size in funzione dei beta
                ax_m_FSS.errorbar(tm, mmag, err_mmag, linestyle = '', marker ='.',label=f'Nlatt = {Nlatt}')


                #stampo su file il max della suscettività che mi serve poi in un altro file per la stima del nu

                dbeta_max =  ( beta[np.argmax(susce)+1] - beta[np.argmax(susce)-1] )/2

                with open(filenu, 'a+' ) as g:
                    print("%d  %f  %f %f %f" % (Nlatt, np.max(susce), err_susce[np.argmax(susce)], beta[np.argmax(susce)], dbeta_max ), file = g)



    ## ANALISI BINDER

                #Funzione di interpolazione per trovare l'intersezione dei binder per Nlatt diversi
                f = interp1d(beta, binder)

                # Grafico del binder bootstrappato in funzione dei beta

                ax_b.errorbar(beta, binder, err_binder, linestyle = '', marker ='.',label=f'Nlatt = {Nlatt}')

                beta, binder = zip(*sorted(zip(beta, binder)))

                ax_b.plot(beta, f(beta))

                # FIGURA ZOOMATA

                # plot the zoomed portion
                betac_pos=beta.index(0.441)

                x_field=beta[betac_pos-1:betac_pos+2]
                y_field=binder[betac_pos-1:betac_pos+2]
                dy_field=err_binder[betac_pos-1:betac_pos+2]

                axins3.errorbar(x_field, y_field, dy_field, marker='.')
                axins3.plot(x_field, f(x_field))



    #ax_c.title.set_text('MODELLO DI ISING 2D \n Calore specifico al variare del beta')
    ax_c.set_xlabel(r'$\beta$')
    ax_c.set_ylabel(r'C')
    ax_c.grid()
    ax_c.minorticks_on()
    ax_c.legend()

    #ax_s.title.set_text('MODELLO DI ISING 2D \n Suscettività al variare del beta')
    ax_s.set_xlabel(r'$\beta$')
    ax_s.set_ylabel(rf'$\chi$')
    ax_s.grid()
    ax_s.minorticks_on()
    ax_s.legend()

    #ax_m.title.set_text('MODELLO DI ISING 2D \n Magnetizzazione intorno alla transizione')
    ax_m.set_xlabel(r'$\beta$')
    ax_m.set_ylabel(r'$\langle |M| \rangle$')
    ax_m.grid()
    ax_m.minorticks_on()
    ax_m.legend()



    ax_b.set_xlabel(r'$\beta$')
    ax_b.set_ylabel(r'$\frac{\langle M^4 \rangle}{\langle M^2 \rangle ^2}$')

    ax_b.grid()
    ax_b.minorticks_on()
    ax_b.legend(loc='lower left')


    ####### GRAFICI per FSS ##########

    #ax_s_FSS.title.set_text('Finite Size Scaling per la Suscettività')
    ax_s_FSS.set_xlabel(r'$\beta-\beta_C$')
    ax_s_FSS.set_ylabel(r'$\chi/L^{\gamma/\nu}$')
    ax_s_FSS.grid()
    ax_s_FSS.minorticks_on()
    ax_s_FSS.legend()

    #ax_m_FSS.title.set_text('Finite Size Scaling per la Magnetizzazione')
    ax_m_FSS.set_xlabel(r'$\beta-\beta_C$')
    ax_m_FSS.set_ylabel(r'$M/L^{-\beta/\nu}$')
    ax_m_FSS.grid()
    ax_m_FSS.minorticks_on()
    ax_m_FSS.legend()



    plt.show()

