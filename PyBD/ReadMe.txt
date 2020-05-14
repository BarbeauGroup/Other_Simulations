------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------
-------------------------- PyBD Python Based Inverse Beta Decay --------------------------
------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------

Last Updated 07/18/2018 - C. Awe, Duke University

PyBD was written using the following version numbers:

ROOT 6.14/00
python 2.7.10
numpy 1.14.3
scipy 1.1.0
matplotlib 2.2.2
pandas 0.22.0
sklearn 0.19.1
tqdm 4.23.4

If python cannot import ROOT, try adding the following lines to your bash_profile (mac):

export PYTHONDIR=/usr/bin/python
export ROOTSYS=/Applications/root_v6.08.06
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib/:$LD_LIBRARY_PATH
export PATH=$ROOTSYS/bin:$PYTHONDIR/bin:$PATH
export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH

Known Issues:
With ROOT 6.08/06 PyRoot cannot write branches properly.

Spectra generated based on:

Reactor:"Calculation of the antineutrino spectrum from thermal fission of 235 U" - Klapdor et. al.
Supernova (Normal Mass Hierarchy):"Neutrino mass hierarchy and three-flavor spectral splits of supernova neutrinos" - Dasgupta et. al. Fig. 2
Supernova (Inverted Mass Hierarchy):"Neutrino mass hierarchy and three-flavor spectral splits of supernova neutrinos" - Dasgupta et. al. Fig. 7

PyRate.py should be kept in the same directory as PyBD.py to work properly.


