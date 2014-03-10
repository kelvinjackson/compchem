http://paton.chem.ox.ac.uk
compchem
========

A collection of scripts useful for computational chemistry...

black_box_RRHO.py
=======
If you are using gaussian to calculate free energies for molecules which are a bit floppy / large then read on…

The main issue here is that the entropy of a harmonic oscillator tends to infinity as the frequency tends to zero - i.e. the harmonic approximation for very low frequencies is a poor one.
This little script that can be applied to the results of gaussian calculations - for vibrational frequencies below a certain (specified) value it replaces the entropy computed for a harmonic oscillator (usually the default in gaussian) by assuming the normal mode acts like a rigid-rotor (http://en.wikipedia.org/wiki/Rigid_rotor) instead. The two regimes are connected using a feathering/damping function so it's not a discontinuous switch between the two. The method is identical to the one used in Grimme, S. Chem. Eur. J. 2012, 18, 9955   

Example of usage:
python black_box_RRHO.py 100 298 PhH_MeCO_FC.out

o  RRHO Free energies using cutoff frequency (cm-1): 100.0  at temperature (K): 298.0
o  PhH_MeCO_FC -385.059930439 0.14874 -384.901665439 -384.944577164 0.583756988422

Shows SCF energy, ZPE, H, G all in Hartree and finally correction (kcal/mol) relative to a fully harmonic G - the corrected value is always higher as previously -TS is always too large.

