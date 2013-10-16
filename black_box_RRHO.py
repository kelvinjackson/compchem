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
AVOGADRO_CONSTANT = 6.0221415e23
AMU_to_KG = 1.66053886E-27
autokcal = 627.509541
kjtokcal = 4.184
atmos = 101.325

# translational energy evaluation (depends on temp.)
def calc_translational_energy(temperature):
	"""
	Calculates the translational energy (kcal/mol) of an ideal gas - i.e. 
	non-interactiing molecules so molar energy = Na * atomic energy
	This approximxation applies to all energies and entropies computed within
	Etrans = 3/2 RT!
	"""
	energy = 1.5 * GAS_CONSTANT * temperature
	energy = energy/kjtokcal/1000.0
	return energy

# rotational energy evaluation (depends on molecular shape and temp.)
def calc_rotational_energy(zpe, symmno, temperature):
	"""
	Calculates the rotaional energy (kcal/mol)
	Etrans = 0 (atomic) ; RT (linear); 3/2 RT (non-linear)
	"""
	if zpe == 0.0: energy = 0.0
	if symmno == 2: energy = GAS_CONSTANT * temperature
	else: energy = 1.5 * GAS_CONSTANT * temperature
	energy = energy/kjtokcal/1000.0
	return energy

# vibrational energy evaluation (depends on frequencies, temp and scaling factor)
def calc_vibrational_energy(frequency_wn, temperature,freq_scale_factor):
	"""
	Calculates the vibrational energy contribution (kcal/mol)
	Includes ZPE (0K) and thermal contributions
	Evib = R * Sum(0.5 hv/k + 1/(e^(hv/KT)-1))
	"""
	energy = 0.0
	frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
	for entry in frequency:
		factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT*temperature))
		temp = factor*temperature*(0.5 + (1/(math.exp(factor)-1)))
		temp = temp*GAS_CONSTANT
		energy = energy + temp
	energy = energy/kjtokcal/1000.0
	return energy

# vibrational Zero point energy evaluation (depends on frequencies and scaling factor)
def calc_zeropoint_energy(frequency_wn,freq_scale_factor):
	"""
	Calculates the vibrational ZPE (kcal/mol)
	EZPE = Sum(0.5 hv/k)
	"""
	energy = 0.0
	frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
	for entry in frequency:
		factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT))
		temp = 0.5*factor
		temp = temp*GAS_CONSTANT
		energy = energy + temp
	energy = energy/kjtokcal/1000.0
	return energy

# translational entropy evaluation (depends on mass, concentration, temp)
def calc_translational_entropy(molecular_mass, conc, temperature):
	"""
	Calculates the translational entropic contribution (cal/(mol*K)) of an ideal gas
	needs the molecular mass
	Convert mass in amu to kg; conc in mol/l to number per m^3
	Strans = R(Ln(2pimkT/h^2)^3/2(1/C)) + 3/2)
	"""
	simass = molecular_mass*AMU_to_KG
	lmda = ((2.0*math.pi*simass*BOLTZMANN_CONSTANT*temperature)**0.5)/PLANCK_CONSTANT
	Ndens = conc*1000*AVOGADRO_CONSTANT
	#Convert molar volume to free volume due to solvent molecule size
	#e.g. chloroform Vol free is 7.5mL
	#e.g. DMF molarity is 12.9 and 77.442 Ang^3 molar volume, so Vfree = 8(10^27/12.9*Na - 
	#4.26)^3 = 3.94 Ang^3 per molecule = 30.6 mL per L 
	# M.H. Abraham, J. Liszi, J. Chem. Soc. Faraday Trans. 74 (1978) 1604.
	
	Ndens = Ndens / (30.6/1000.0)
	entropy = GAS_CONSTANT*(2.5+math.log(lmda**3/Ndens))/4.184
	
	return entropy


# rotational entropy evaluation (depends on molecular shape and temp.)
def calc_rotational_entropy(zpe, symmno, rotemp, temperature):
	"""
		Calculates the rotational entropy (cal/(mol*K))
		Strans = 0 (atomic) ; R(Ln(q)+1) (linear); R(Ln(q)+3/2) (non-linear)
		"""
	#print moment_of_inertia

	#rotemp = [PLANCK_CONSTANT**2/(8*math.pi**2*BOLTZMANN_CONSTANT*entry*AMU_to_KG*1E-20) for entry in moment_of_inertia]
	#print rotemp
	qrot = math.pi*temperature**3/(rotemp[0]*rotemp[1]*rotemp[2])
	qrot = qrot ** 0.5
	qrot = qrot/symmno
	if zpe == 0.0: entropy = 0.0
	if symmno == 2: entropy = GAS_CONSTANT * (math.log(qrot/2.0) + 1)
	else: energy = entropy = GAS_CONSTANT * (math.log(qrot) + 1.5)
	entropy = entropy/kjtokcal
	return entropy


# harmonic entropy evaluation
def calc_harmonic_entropy(frequency_wn, temperature,freq_scale_factor):
	"""
	Calculates the entropic contribution (cal/(mol*K)) of a harmonic oscillator for
	a list of frequencies of vibrational modes
	Sv = R(hv/(k(e^(hv/KT)-1) - ln(1-e^(-hv/kT)))	
	"""
	entropy = []
	frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
	for entry in frequency:
		factor = ((PLANCK_CONSTANT*entry*freq_scale_factor)/(BOLTZMANN_CONSTANT*temperature))
		temp = factor*(1/(math.exp(factor)-1)) - math.log(1-math.exp(-factor))
		temp = temp*GAS_CONSTANT/4.184
		entropy.append(temp)
	return entropy


# rigid  rotor entropy evaluation
def calc_rrho_entropy(frequency_wn, temperature,freq_scale_factor):
	"""
	Calculates the entropic contribution (cal/(mol*K)) of a rigid-rotor harmonic oscillator for
	a list of frequencies of vibrational modes
	Sr = R(1/2 + 1/2ln((8pi^3u'kT/h^2))
	"""
	Bav = 10.0e-44
	entropy = []
	frequency = [entry * SPEED_OF_LIGHT for entry in frequency_wn]
	
	for entry in frequency:
		mu = PLANCK_CONSTANT/(8*math.pi**2*entry*freq_scale_factor)
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
			if line.strip().startswith('Molecular mass:'): molecular_mass = float(line.strip().split()[2])
			if line.strip().startswith('Rotational symmetry number'): symmno = int((line.strip().split()[3]).split(".")[0])
			if line.strip().startswith('Rotational temperatures'): rotemp = [float(line.strip().split()[3]), float(line.strip().split()[4]), float(line.strip().split()[5])]
		 
				
		# Calculate Translational, Rotational and Vibrational contributions to the energy
		Utrans = calc_translational_energy(temperature)
		Urot = calc_rotational_energy(self.zero_point_corr, symmno, temperature)
		Uvib = calc_vibrational_energy(frequency_wn, temperature,freq_scale_factor)
		ZPE = calc_zeropoint_energy(frequency_wn, freq_scale_factor)
		
		# Calculate Translational, Rotational and Vibrational contributions to the entropy
		Strans1atm = calc_translational_entropy(molecular_mass, atmos/(GAS_CONSTANT*temperature), temperature)
		Strans = calc_translational_entropy(molecular_mass, conc, temperature)
		conc_correction = Strans - Strans1atm
			
		Srot = calc_rotational_entropy(self.zero_point_corr, symmno, rotemp, temperature)
				
		# Calculate harmonic entropy, free-rotor entropy and damping function for each frequency - functions defined above
		Svibh = calc_harmonic_entropy(frequency_wn, temperature,freq_scale_factor)
		Svibrrho = calc_rrho_entropy(frequency_wn, temperature,freq_scale_factor)
		damp = calc_damp(frequency_wn, FREQ_CUTOFF)
		
		# Compute entropy (cal/mol/K) using the two values and damping function
		vib_entropy = []
		for j in range(0,len(frequency_wn)):
			vib_entropy.append(Svibh[j] * damp[j] + (1-damp[j]) * Svibrrho[j])

		Svib = sum(vib_entropy)
		RRHO_correction = Svib - sum(Svibh)
		
		# Add all terms to get Free energy

		self.enthalpy = self.scf_energy + (Utrans + Urot + Uvib + GAS_CONSTANT*temperature/kjtokcal/1000.0)/autokcal
		self.zpe = ZPE/autokcal
		self.entropy = (Strans + Srot + Svib)/autokcal/1000.0
		self.gibbs_free_energy = self.enthalpy - temperature * self.entropy

		self.RRHO_correction = -RRHO_correction * temperature/1000.0
		self.conc_correction = -conc_correction * temperature/1000.0

if __name__ == "__main__":
	
	# Takes arguments: cutoff_freq g09_output_files
	files = []
	if len(sys.argv) > 3:
		FREQ_CUTOFF = float(sys.argv[1])
		temperature = float(sys.argv[2])
		conc = float(sys.argv[3])
		freq_scale_factor = float(sys.argv[4])
		
		for i in range(5,len(sys.argv)):
			files.append(sys.argv[i])
	else:
		print "\nWrong number of arguments used. Correct format: black_box_entropy cutoff_freq temp concn g09_output_files\n"
		sys.exit()
	
	print "\no  RRHO Free energies using cutoff frequency (cm-1):",FREQ_CUTOFF, " at temperature (K):", temperature, " and conc. (mol/l)", conc
	for file in files:
		bbe = calc_bbe(file, FREQ_CUTOFF, temperature)
		print "o ",file.split(".")[0],
		if not hasattr(bbe,"gibbs_free_energy"): print "Warning! Job did not finish normally!"
		if hasattr(bbe, "scf_energy"): print bbe.scf_energy,
		else: print "N/A",
		if hasattr(bbe, "zero_point_corr"): print bbe.zpe,
		else: print "N/A",
		if hasattr(bbe, "enthalpy"): print bbe.enthalpy,
		else: print "N/A",
		if hasattr(bbe, "gibbs_free_energy"): print bbe.gibbs_free_energy,		
		if hasattr(bbe, "RRHO_correction"): print bbe.RRHO_correction,
		else: print ""
		if hasattr(bbe, "conc_correction"): print bbe.conc_correction
		
