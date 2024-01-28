Study of Ising Model 2D numerically.

To get the simulated data run simulatore.c (needed simulazione.h, input.txt and beta.txt). 

BOOTSTRAP_ISING.c makes a binned bootstrap of the sample created previously to get the proper error on the quantities analyzed.

bootstrap_per_ck.c is used to make bootstrap for data created by CK_CALCOLO.py aimed to verify the correalations in data genrated by the algorithm and plotted
by C_kappa_view.py.

Other Python files are used for the analysis and plot of the interesting observables. Everything is fully explained in .pdf
and inside files themselves.
