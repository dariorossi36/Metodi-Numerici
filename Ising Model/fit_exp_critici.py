# cartella principale
root_dir = "/home/dario/Documents/UNI/Metodi/Modulo1/Ising/"

import pylab
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit


# Dati reticolo
Nlatt = 40


# Definizione del beta critico (valore numerico da trovare con plot di binder in funzione di L)
beta_C = 0.441
beta_c_err = 0.005


# Importo dati del bootstrap
filename = root_dir + f"Nlatt={Nlatt}" + "/Bootstrap/bootstrap.txt"

ene, err_ene, mag, err_mag, calore, err_calore, susce, err_susce, binder, err_binder, beta = np.loadtxt(filename, unpack=True)

print(f"\n \n Nlatt={Nlatt} \n \n")

################################################################################################
#################### FIT ESPONENTE CRITICO alpha (CALORE SPECIFICO vs t) #######################
################################################################################################

# Tagliamo le variabili per fare il fit (quindi prendendo da quando scende la curva)
inizio1 = int(np.argwhere(beta==beta_C))

t = (np.abs(beta - beta_C))    # temperatura ridotta

inizio1 = inizio1 - 17
if(Nlatt>=50):
    inizio1 = inizio1 - 38

x=t[:inizio1]
y1=calore[:inizio1]
dy1=err_calore[:inizio1]
dx = np.ones(len(x))*beta_c_err

# Rimozione di altri eventuali outlier
# outlier=8
# x=np.delete(x,outlier)
# y1=np.delete(y1,outlier)
# dy1=np.delete(dy1,outlier)
# dx=np.delete(dx,outlier)

# Funzione di fit (uguale per tutti, si defisce solo qua)
def f(x,a,b,c):
    return a*x**(b)+c

# Procedura di best fit
init=(1,0,0)
popt, pcov=curve_fit(f, x, y1, init, dy1, maxfev = 100000000)

ndof=len(x)-len(init)
chi2=(((y1-f(x, *popt))/dy1)**2).sum()
print('************** CALORE SPECIFICO *************')
print('passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=dy1
for i in range(0, 3, 1):
    dxy=np.sqrt(dxy**2+dx**2)
    popt, pcov=curve_fit(f, x, y1, popt, dxy,maxfev = 100000000)
    chi2=(((y1-f(x, *popt))/dxy)**2).sum()
    print('passo %d' % i)
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))


print('\n\n\n')


# Grafico della funzione di fit

xx=np.linspace(min(x), max(x)+0.0005, 1000)

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,8), height_ratios=[2,1])

ax1.set_title(f'Nlatt={Nlatt}', size=13, fontweight='bold')
ax1.errorbar(x, y1, dy1,linestyle = '', marker = '.')
ax1.plot( xx, f(xx, *popt))
ax1.set_ylabel('C')

ax1.grid(color='lightgray')
ax1.minorticks_on()


#residui normalizzati

r = (y1-f(x,*popt))/dxy
ax2.errorbar( x,r, linestyle='', marker='.')

ax2.set_xlabel(r'$|\beta-\beta_C|$')
ax2.set_ylabel('residui norm')
ax2.grid(color='lightgray')

# salvo valori best fit in un file
#print("%f  %f  %f" % (Nlatt, popt[1], np.sqrt(pcov.diagonal())[1]), file = f)

################################################################################################
####################### FIT ESPONENTE CRITICO beta (SUSCETTIVITÀ vs t) ########################
################################################################################################

# Tagliamo le variabili per fare il fit (quindi prendendo da quando scende la curva)
inizio = int(np.argwhere(beta==beta_C))
inizio = inizio + 8

if(Nlatt>40):
    inizio = inizio -3

x=t[inizio:]
y=susce[inizio:]
dy=err_susce[inizio:]
dx = np.ones(len(x))*beta_c_err

outlier=(0,1,2)
x=np.delete(x,outlier)
y=np.delete(y,outlier)
dy=np.delete(dy,outlier)
dx=np.delete(dx,outlier)



# Procedura di best fit
init=(1,-7/4,0)
popt, pcov=curve_fit(f, x, y, init, dy,maxfev = 100000000)


ndof=len(x)-len(init)
chi2=(((y-f(x, *popt))/dy)**2).sum()
chi2=(((y-f(x, *popt))/dy)**2).sum()
print('************** SUSCETTIVITA *************')
print('passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=dy
for i in range(0, 3, 1):
    dxy=np.sqrt(dxy**2+dx**2)
    popt, pcov=curve_fit(f, x, y, popt, dxy, maxfev = 100000000)
    chi2=(((y-f(x, *popt))/dxy)**2).sum()
    print('passo %d' % i)
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))


print('\n\n\n')


# Grafico della funzione di fit

xx=np.linspace(min(x), max(x)+0.0005, 1000)
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,8), height_ratios=[2,1])

ax1.set_title(f'Nlatt={Nlatt}', size=13, fontweight='bold')
ax1.errorbar(x, y, dy,linestyle = '', marker = '.')
ax1.plot( xx, f(xx, *popt))
ax1.set_ylabel(r'$\chi$')

ax1.grid(color='lightgray')
ax1.minorticks_on()


#residui normalizzati

r = (y-f(x,*popt))/dxy
ax2.errorbar( x,r, linestyle='', marker='.')
ax2.set_xlabel(r'$ |t| = \beta-\beta_C$')
ax2.set_ylabel('residui norm')
ax2.grid(color='lightgray')

# salvo valori best fit in un file
#print("%f  %f" % (popt[1], np.sqrt(pcov.diagonal())[1]), file = f)

################################################################################################
############## FIT ESPONENTE CRITICO gamma (MAGNETIZZAZIONE vs t) ##############################
################################################################################################

# Tagliamo le variabili per fare il fit (quindi prendendo da quando scende la curva)
#inizio2 = np.argmax(magne)

t = (np.abs(beta - beta_C))    # temperatura ridotta

inizio2=int(np.argwhere(beta==beta_C)) +4

if(Nlatt<20):
    inizio2 = inizio2 + 4

x=t[inizio2:]
y2=mag[inizio2:]
dy2=err_mag[inizio2:]
dx = np.ones(len(x))*beta_c_err

# outlier= len(x)-1,len(x)-2
# x=np.delete(x,outlier)
# y2=np.delete(y2,outlier)
# dy2=np.delete(dy2,outlier)
# dx=np.delete(dx,outlier)


init=(-1,1/8,0)
popt, pcov=curve_fit(f, x, y2, init, dy2, maxfev = 100000000)
# m_fit, q_fit=popt
# dm_fit, dq_fit=np.sqrt(pcov.diagonal())
# print('primo fit m=%f +- %f' %(m_fit, dm_fit ) )


ndof=len(x)-len(init)
chi2=(((y2-f(x, *popt))/dy2)**2).sum()
chi2=(((y2-f(x, *popt))/dy2)**2).sum()
print('************** MAGNETIZZAZIONE *************')
print('passo zero')
print('popt:', popt)
print('dpopt:', np.sqrt(pcov.diagonal()))
print('chi2=%f, \nndof=%f' %(chi2, ndof))
dxy=dy2
for i in range(0, 5, 1):
    dxy=np.sqrt(dxy**2+dx**2)
    popt, pcov=curve_fit(f, x, y2, popt, dxy, maxfev = 100000000)
    chi2=(((y2-f(x, *popt))/dxy)**2).sum()
    print('passo %d' % i)
    print('popt:', popt)
    print('dpopt:', np.sqrt(pcov.diagonal()))
    print('chi2=%f, \nndof=%f' %(chi2, ndof))


print('\n\n\n')


# Grafico della funzione di fit

xx=np.linspace(min(x), max(x)+0.0005, 1000)
fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,8), height_ratios=[2,1])

ax1.set_title(f'Nlatt={Nlatt}', size=13, fontweight='bold')
ax1.errorbar(x, y2, dy2,linestyle = '', marker = '.')
ax1.plot( xx, f(xx, *popt))
#ax1.xscale('log')
ax1.set_ylabel(rf'$\langle |M| \rangle$')
ax1.grid(color='lightgray')
ax1.minorticks_on()


#residui normalizzati

r = (y2-f(x,*popt))/dxy
ax2.errorbar( x,r, linestyle='', marker='.')
ax2.set_ylabel('residui norm')
ax2.set_xlabel(r'$ |t| = \beta-\beta_C$')
ax2.grid(color='lightgray')

# salvo valori best fit in un file
#print("%f  %f \n" % (popt[1], np.sqrt(pcov.diagonal())[1] ), file = f)

# ###############################################################################################
# #################### FIT ESPONENTE CRITICO nu (CHI_max vs L) ##################################
# ###############################################################################################

# x,y,dy,beta_max, dbeta_max=np.loadtxt(root_dir+"/Exp_crit/exp_nu.txt", unpack=True)

# ## FIT PER TROVARE gamma/nu
# def f(x,a,b,c):
#     return a*x**(b)+c


# init=(1,-1.75,0)

# popt, pcov=curve_fit(f, x, y, init, dy)
# ndof = len(x)-len(init)
# chi2=(((y-f(x, *popt))/dy)**2).sum()
# print('************** SUSCETTIVTÀ *************')
# print('passo 0')
# print('popt:', popt)
# print('dpopt:', np.sqrt(pcov.diagonal()))
# print('chi2=%f, \nndof=%f' %(chi2, ndof))
# dx=0
# dxy=dy
# for i in range(1, 4, 1):
#     dxy=np.sqrt(dxy**2+dx**2)
#     popt, pcov=curve_fit(f, x, y, popt, dxy)
#     chi2=(((y-f(x, *popt))/dxy)**2).sum()
#     print('passo %d' % i)
#     print('popt:', popt)
#     print('dpopt:', np.sqrt(pcov.diagonal()))
#     print('chi2=%f, \nndof=%f' %(chi2, ndof))


# print('\n\n\n')

# #figura

# xx=np.linspace(min(x), max(x)+0.0005, 1000)
# fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,8), height_ratios=[2,1])

# ax1.errorbar(x, y, dy,linestyle = '', marker = '.')
# ax1.plot( xx, f(xx, *popt))

# ax1.set_ylabel(r'$\chi_{max}$')
# #plt.title('MODELLO DI ISING 2D \n Suscettività massima in funzione di N')
# ax1.grid(color='lightgray')
# ax1.minorticks_on()

# #residui normalizzati

# r = (y-f(x,*popt))/dxy
# ax2.errorbar( x,r, linestyle='', marker='.')

# ax2.set_xlabel(r'$N_{latt}$')
# ax2.set_ylabel('residui norm')
# ax2.grid(color='lightgray')

# ## FIT PER TROVARE ANCHE IL beta_critico

# def f(x,x_bar,nu,beta_c):
#     return x_bar*x**(-1/nu)+beta_c


# init=(1,1,0.441)
# # cambio variabile per comodità
# y=beta_max
# dy = dbeta_max
# popt, pcov=curve_fit(f, x, y, init, dy)
# ndof = len(x)-len(init)
# chi2=(((y-f(x, *popt))/dy)**2).sum()
# print('************** SUSCETTIVTÀ + beta_crit *************')
# print('passo 0')
# print('popt:', popt)
# print('dpopt:', np.sqrt(pcov.diagonal()))
# print('chi2=%f, \nndof=%f' %(chi2, ndof))
# dx=0
# dxy=dy
# for i in range(1, 4, 1):
#     dxy=np.sqrt(dxy**2+dx**2)
#     popt, pcov=curve_fit(f, x, y, popt, dxy)
#     chi2=(((y-f(x, *popt))/dxy)**2).sum()
#     print('passo %d' % i)
#     print('popt:', popt)
#     print('dpopt:', np.sqrt(pcov.diagonal()))
#     print('chi2=%f, \nndof=%f' %(chi2, ndof))


# print('\n\n\n')

# #figura

# xx=np.linspace(min(x), max(x)+0.0005, 1000)
# fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,8), height_ratios=[2,1])

# ax1.errorbar(x, y, dy,linestyle = '', marker = '.')
# ax1.plot( xx, f(xx, *popt))

# ax1.set_ylabel(r'$\beta_{pc}$')
# #plt.title('MODELLO DI ISING 2D \n Suscettività massima in funzione di N')
# ax1.grid(color='lightgray')
# ax1.minorticks_on()

# #residui normalizzati

# r = (y-f(x,*popt))/dxy
# ax2.errorbar( x,r, linestyle='', marker='.')

# ax2.set_xlabel(r'$N_{latt}$')
# ax2.set_ylabel('residui norm')
# ax2.grid(color='lightgray')

plt.show()