#!/usr/bin/python

###                   ###                     ###          ###      ###
###                   ###                     ###          ###      ###
    #####b.   ####b.  ###### .d##b.  #####b.  ###  ####b.  #####b.  ###
### ### "##b     "##b ###   d##""##b ### "##b ###     "##b ### "##b ###
### ###  ### .d###### ###   ###  ### ###  ### ### .d###### ###  ###
### ### d##P ###  ### Y##b. Y##..##P ###  ### ### ###  ### ### d##P ###
### #####P"  "Y######  "Y### "Y##P"  ###  ### ### "Y###### #####P"  ###
    ###                                                             
    ###

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk

#######################################################################
#                      black_box_RRHO.py                              #
#  A program to computes the vibrational entropy using                #
#  a harmonic approximation for frequencies above a certain           #
#  cutoff, and a free-rotor approximation below. A damping            #
#  function feathers the entropy term between the two limits          #
#  as described in: Grimme, S. Chem. Eur. J. 2012, 18, 9955           #
#######################################################################
#######  Written by:  Rob Paton #######################################
#######  Last modified:  Mar 20, 2013 #################################
#######################################################################

## Adapted from Dr Arne Dieckmann's (UCLA) quasiharmonic free energy python code##

import sys
import math

# CONSTANTS
GAS_CONSTANT = 8.3144621
PLANCK_CONSTANT = 6.62606957e-34 
BOLTZMANN_CONSTANT = 1.3806488e-23 
SPEED_OF_LIGHT = 2.99792458e10

# harmonic entropy evaluation 
def calc_harmonic_entropy(frequency_wn, temperature):
	"""
	Calculates the entropic contribution (cal/(mol*K)) of a harmonic oscillator for
	a list of frequencies of vibrational modes
	Sv = R(hv/(k(e^(hv/KT)-1) - ln(1-e^(-hv/kT)))	
	"""
	entropy = []
	frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
	for entry in frequency:
		factor = ((PLANCK_CONSTANT*entry)/(BOLTZMANN_CONSTANT*temperature))
		temp = factor*(1/(math.exp(factor)-1)) - math.log(1-math.exp(-factor))
		temp = temp*GAS_CONSTANT/4.184
		entropy.append(temp)
	return entropy


# rigid  rotor entropy evaluation
def calc_rrho_entropy(frequency_wn, temperature):
	"""
	Calculates the entropic contribution (cal/(mol*K)) of a rigid-rotor harmonic oscillator for
	a list of frequencies of vibrational modes
	Sr = R(1/2 + 1/2ln((8pi^3u'kT/h^2))
	"""
	Bav = 10.0e-44
	entropy = []
	frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
	
	for entry in frequency:
		mu = PLANCK_CONSTANT/(8*math.pi**2*entry)
		muprime = mu*Bav/(mu +Bav)
		#print entry, mu,muprime
		
		factor = (8*math.pi**3*muprime*BOLTZMANN_CONSTANT*temperature)/(PLANCK_CONSTANT**2)
		temp = 0.5 + math.log(factor**0.5)
		temp = temp*GAS_CONSTANT/4.184
		#print temp
		entropy.append(temp)
	return entropy


# damping function
def calc_damp(frequency_wn, FREQ_CUTOFF):
	"""
	Calculates the Head-Gordon damping function with alpha=4
	"""
	alpha = 4
	damp = []
	for entry in frequency_wn:
		omega = 1/(1+(FREQ_CUTOFF/entry)**alpha)
		damp.append(omega)
	return damp


class calc_bbe:	
	def __init__(self, file, FREQ_CUTOFF, temperature):

		# Frequencies in waveunmbers
		frequency_wn = []
		
		# Read commandline arguments
		g09_output = open(file, 'r')
		
		# Iterate over output
		for line in g09_output:
			# look for low frequencies  
			if line.strip().startswith('Frequencies --'):
				for i in range(2,5):
					x = float(line.strip().split()[i])
					#  only deal with real frequencies
					if x > 0.00: frequency_wn.append(x)
			# look for SCF energies, last one will be correct
			if line.strip().startswith('SCF Done:'): self.scf_energy = float(line.strip().split()[4])
			# look for thermal corrections 
			if line.strip().startswith('Zero-point correction='): self.zero_point_corr = float(line.strip().split()[2])
			if line.strip().startswith('Thermal correction to Energy='): self.energy_corr = float(line.strip().split()[4])
			if line.strip().startswith('Thermal correction to Enthalpy='): enthalpy_corr = float(line.strip().split()[4])
			if line.strip().startswith('Thermal correction to Gibbs Free Energy='): gibbs_corr = float(line.strip().split()[6])

		# Calculate harmonic entropy, free-rotor entropy and damping function for each frequency - functions defined above
		Sv = calc_harmonic_entropy(frequency_wn, temperature)
		Sr = calc_rrho_entropy(frequency_wn, temperature)
		damp = calc_damp(frequency_wn, FREQ_CUTOFF)
		
		# Compute entropy (cal/mol/K) using the two values and damping function
		vib_entropy = []
		for j in range(0,len(frequency_wn)):
			vib_entropy.append(Sv[j] * damp[j] + (1-damp[j]) * Sr[j])
				
		# The difference between this adjusted Svib and the purely harmonic one is used to correct the Free energy
		# Converted to kcal/mol
		self.TSdiff = (sum(Sv)-sum(vib_entropy))*temperature/1000.0
		#print self.TSdiff
	
		self.new_gibbs_corr = gibbs_corr + self.TSdiff/627.51
		self.enthalpy = self.scf_energy + enthalpy_corr
		self.harmonic_gibbs_energy = self.scf_energy + gibbs_corr
		self.gibbs_energy = self.scf_energy + self.new_gibbs_corr
		self.entropy = (self.enthalpy - self.gibbs_energy)/temperature
		

if __name__ == "__main__":
	
	# Takes arguments: cutoff_freq g09_output_files
	files = []
	if len(sys.argv) > 3:
		FREQ_CUTOFF = float(sys.argv[1])
		temperature = float(sys.argv[2])
		for i in range(3,len(sys.argv)):
			files.append(sys.argv[i])
	else:
		print "\nWrong number of arguments used. Correct format: black_box_entropy cutoff_freq temp g09_output_files\n"
		sys.exit()
	
	print "\no  RRHO Free energies using cutoff frequency (cm-1):",FREQ_CUTOFF, " at temperature (K):", temperature
	for file in files:
		bbe = calc_bbe(file, FREQ_CUTOFF, temperature)
		print "o ",file.split(".")[0],
		if not hasattr(bbe,"gibbs_energy"): print "Warning! Job did not finish normally!"
		if hasattr(bbe, "scf_energy"): print bbe.scf_energy,
		else: print "N/A",
		if hasattr(bbe, "zero_point_corr"): print bbe.zero_point_corr,
		else: print "N/A",
		if hasattr(bbe, "enthalpy"): print bbe.enthalpy,
		else: print "N/A",
		if hasattr(bbe, "gibbs_energy"): print bbe.gibbs_energy, bbe.TSdiff
		else: print "N/A"
