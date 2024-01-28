Files to study a scalar field in 2 dimension with a path integral approach.

To get the simulated data run scalar2d.c (needed scalar2d.h). It will be asked to choose betwwn two different simulations: "limite al continuo",
which uses valori_cont.txt and make the simulation at a fixed temperature, varying the lattice spacing towards continuity; "Temperatura" (uses valori_temp.txt),
that make the simulation at a fixed value of lattice spacing but varying temperature.

BOOTSTRAP.c makes a binned bootstrap of the sample created previously to get the proper error on the quantities analyzed.

Python files are used for the analysis and plot of the interesting observables. everything is fully explained in .pdf and inside files themselves.
