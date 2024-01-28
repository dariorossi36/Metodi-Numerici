Study of a quantum harmonic oscillator with a path integral approach.

To get the simulated data run oscillazione_eta.c (needed oscillazione.h and input.txt). It will be asked to choose between two different simulations: 
"limite al continuo", which uses valori_netafisso.txt and make the simulation at a fixed temperature, varying the lattice spacing towards continuity; 
"Temperatura" (uses valori_n.txt), that make the simulation at a fixed value of lattice spacing but varying temperature.

BOOTSTRAP_OSC.c makes a binned bootstrap of the sample created previously to get the proper error on the quantities analyzed.

oscillazione_eta_splitting.c is used to make simulations for the staudy of the splitting of energy levels (uses splitting_energy.h and valori_split.txt)
and also for the plot of the fundamental state function (funz_fondamentale.py). To see the splitting the order of run is: oscillazione_eta_splitting.c ->
 CK_CALCOLO2.py -> CK_BOOTSTRAP.c -> splitting_view.py.

Other Python files are used for the analysis and plot of the interesting observables. Everything is fully explained in Relazione_oscillatore.pdf
and inside files themselves.
