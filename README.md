# sigma-trick-reference
Reference implementations of the sigma trick for various symmetry cases.

Each of the folders corresponds to a different symmetry type as indicated by the file name.

The implementation of the self-energy trick is in the file sigmatrick-new.jl

# Usage

Each folder contains the spectral data in the subfolder spectral_Nz# -- this has already been broadened from raw
spectral NRG data. If you want to redo the NRG calculation using [NRG Ljubljana](http://nrgljubljana.ijs.si), the 
appropriate param.loop is provided in every folder, as well as 1_run script. The data are generated by running this
script. Note that you also need a working Mathematica installation with the [SNEG package](http://nrgljubljana.ijs.si/sneg/) 
installed.

The calculated correlation functions and self-energies are in the files: c-Fl.dat, c-Fr.dat, c-G.dat, c-I.dat, aw.dat, aw-old.dat.
The files c-Fl.dat, c-Fr.dat, c-G.dat and c-I.dat are the correlators defined in [Kugler's paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.105.245132),
while aw.dat and aw-old.dat are the imaginary parts of the impurity spectral function, the first calculated with the new self-energies 
(contained in the file sigma.dat) and the latter with the old self-energies (in the file sigma-old.dat).