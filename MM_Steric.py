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
# AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk

#######################################################################
#######################################################################
#######  Written by:  Rob Paton #######################################
#######  Last modified:  Feb 07, 2014 #################################
#######################################################################

interactivemode = 1
E_printcutoff = 2.0
dashedline = "   - - - - - - - - - - - - - - - - - - - - - - - - - - - - "

# For reading Gaussian formatted input/output files
from ccParse import *

#Python libararies
import random, sys, os, commands, string, math, numpy

## Check for integer when parsing ##
def is_number(s):
    try: int(s); return True
    except ValueError: return False

## Arrays for attractive and repulsive interactions ##
attractive_vdw=[0]
repulsive_vdw=[0]
total_vdw=[0]

## Conversion factors ##
autoang = 0.52917726
autokcal = 627.509541
Eps = 332.0522173

## Work out the atoms in the same molecules
def getMollist(bondmatrix,startatom,natom):
	
	# The list of atoms in a molecule
	atomlist=[]
	atomlist.append(startatom)
	molecule1=[]
	nextlot=[]
	count = 0

	while count<1000:
		nextlot=[]
		#print count,
		for atom in atomlist:
			#print atom, bondmatrix[atom]
			for i in range(0,len(bondmatrix[atom])):
				if bondmatrix[atom][i] == 1:
					alreadyfound = 0
					for at in atomlist:
						if i == at: alreadyfound = 1
					if alreadyfound == 0 and i < natom: atomlist.append(i)

		count=count+1	
	return atomlist
	

## The periodic table...
elements=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

## For a specified elemental symbol return the atomic number
def atomicnumber(element):
	atomicno = 0
	for i in range(0,len(periodictable)):
		if element == periodictable[i]: atomicno = i
	return atomicno


## Convert from G09 format atom types number codes to standard string
## put these into uff.parm
def at_to_string(type):
	numcodes = [0, 10011000, 10151004, 10061000, 10061003, 10081000, 10081002, 10081003]
	at = ["null", "H-H_", "P-P_3+5", "O-O_R", "O-O_3", "C-C_R", "C-C_2", "C-C_3"]
	valid = 0
	for i in range(0,len(numcodes)):
		if numcodes[i] == type: valid = i
	if valid == 0: print type, "not found"
	return at[valid]


## Bondi atomic radius
def bondiRadius(Z):
	#Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452, except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391
	#Radii that are not available in either of these publications have RvdW = 2.00 Angstrom
	bondi = [0.0,1.09, 1.40, 1.82,2.00,2.00,1.70,1.55,1.52,1.47,1.54,2.27,1.73,2.00,2.10,1.80,1.80,1.75,1.88,2.75,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,1.87,2.00,1.85,1.90,
			 1.85,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,2.00,2.06,1.98,2.16,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.72,1.66,1.55,1.96,2.02,2.00,2.00,2.00,
			 2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.86]
	if Z<len(bondi): radius = bondi[Z]
	else: radius = 2.0
	return radius


## UFF parameters for VdW interactions
def getPARM(ff,type):
	parmfile = ff+".prm"
	if not os.path.exists(parmfile): print ("\nFATAL ERROR: "+ff+" parameter file [ %s ] does not exist"%parmfile)
	else:
		parmlines = open(parmfile,"r")
		parms = parmlines.readlines()
	#print parms
	
	for par in parms:
		#print par
		if par.find(type) > -1:
			#print par.split()
			vdw = [float(par.split()[2]),float(par.split()[3])]
			name = par.split()[1]
	try:
		#print vdw
		return vdw
	except NameError:
		print "can't find param!"
		return [0.0,0.0,0.0]


## The computation of the D3 dispersion correction
class calcEDA:
	def __init__(self, comfile):
		
		## Reading Gaussian09 input - connectivity information
		if comfile != "none":
			print "\no  SCRAPING G09 COM FILE <", comfile,">"
			comData = getinData(comfile)
			if hasattr(comData, "CARTESIANS"): print "   Read", len(comData.CARTESIANS), "Atomic Cartesian coordinates",; cartesians = comData.CARTESIANS
			if hasattr(comData, "ATOMTYPES"): print "   Read", len(comData.ATOMTYPES), "atomtypes",; atomtypes = comData.ATOMTYPES; natom = len(atomtypes)
			if hasattr(comData, "CONNECTIVITY"): print "   Read", len(comData.CONNECTIVITY), "Atomic Connectivity Information",; bondedlist = comData.CONNECTIVITY; #connectivity = comData.BONDINDEX
			#if len(comData.ATOMTYPES) != len(comData.LEVELTYPES): print "o  WARNING: problem assigning atom types /ONIOM levels!"; sys.exit(0)
			#if len(comData.ATOMTYPES) != len(comData.MULLIKEN): print "o  WARNING: PDB and G09 file have different numbers of atoms!"; sys.exit(0)
			
			try: bondedlist
			except NameError:
				comData.CONNECTIVITY = []; #comData.BONDINDEX = []
				for i in range(0,len(comData.CARTESIANS)):
					#print i,
					comData.CONNECTIVITY.append([]); #comData.BONDINDEX.append([])
					for j in range(0,len(comData.CARTESIANS)):
						#comData.BONDINDEX[i].append(0)
						if numpy.linalg.norm(numpy.array(comData.CARTESIANS[i])-numpy.array(comData.CARTESIANS[j])) < 0.5*(bondiRadius(atomicnumber(comData.ATOMTYPES[i])) + bondiRadius(atomicnumber(comData.ATOMTYPES[j]))):
							if i != j: comData.CONNECTIVITY[i].append(str(j+1)+"__1.0"); #comData.BONDINDEX[i][j] = 1

				print "\n   Generated atomic connectivity information from Cartesian coordinates"; bondedlist = comData.CONNECTIVITY; #connectivity = comData.BONDINDEX


			try: atomtypes
			except NameError:
				atomtypes = []
				for i in range(0,len(comData.ATOMICTYPES)):
					if comData.ATOMICTYPES[i] == 0 : atomtypes.append(comData.ATOMTYPES[i])
					else: atomtypes.append(at_to_string(comData.ATOMICTYPES[i]))
		
		
		nmols = 1
		print "\n   There are", nmols, "molecules present:"
		#for i in range(0,len(mols)): print "   Molecule", (i+1), "containing",len(mols[i]), "atoms"
		
		# Array for interaction energies by residue
		Electrostatic_term = []
		VdW_term = []


		## Interatomic EDA
		eda_model = 1
		## Compute the pairwise interaction energies between atoms i and j
		if eda_model == 1:
			print "\n                          ", ("Charge").rjust(10), ("E(el)").rjust(10), ("E(vdw)").rjust(10)
			for i in range(0,natom):
				VdW_term.append(0.0); Electrostatic_term.append(0.0)
				for j in range(0,natom):
					if i != j:
						dist = numpy.linalg.norm(numpy.array(comData.CARTESIANS[i])-numpy.array(comData.CARTESIANS[j]))
						
						# If atoms are 1,2- 1,3- or 1,4-separated, scaling factors are used
						Esfactor=1.0; VDWfactor =1.0
						
						#if connectivity[i][j]==1: Esfactor = 0; VDWfactor = 0 #When I and J are directly bonded
						#else:
						for bondediatom in bondedlist[i]:
								for bondedjatom in bondedlist[j]:
									if int(bondediatom.split("__")[0]) == int(bondedjatom.split("__")[0]):
										Esfactor = 0; VDWfactor = 0 # When I and J are bonded to a connecting atom
									
									for bondedkatom in bondedlist[int(bondediatom.split("__")[0])-1]:
										if int(bondedkatom.split("__")[0]) == int(bondedjatom.split("__")[0]):
											if int(bondedkatom.split("__")[0]) != int(bondediatom.split("__")[0]):
												Esfactor = 1/1.2; VDWfactor = 1/2.0 # When I and J are bonded to 2 connecting atoms
						
						# Calculate scaled Electrostatic interaction (factor of 0.5 is because each interaction is counted twice)
						Ees = 0.0 #0.5 * Esfactor*Eps*charges[i]*charges[j]/dist
						#print str(i+1).rjust(6), str(j+1).rjust(6), str(Esfactor).rjust(6), str(charges[i]).rjust(10), str(charges[j]).rjust(10), str(Ees).rjust(10)
						Electrostatic_term[i] = Electrostatic_term[i] + Ees

						mm_type =2
						if mm_type == 1:
							# Calculate  AMBER 6-12 VDW parameters from atom types
							parmi = getPARM("amber", atomtypes[i])
							parmj = getPARM("amber", atomtypes[j])
							Eijstar = (parmi[1]*parmj[1])**0.5
							Rijstar = parmi[0] + parmj[0]
							Rijstar_6 = Rijstar ** 6
							Rijstar_12 = Rijstar_6 ** 2
							Aij = Eijstar * Rijstar_12
							Bij = 2.0 * Eijstar * Rijstar_6

						if mm_type == 2:
							# Calculate UFF 6-12 VDW parameters from atom types
							parmi = getPARM("uff", atomtypes[i])
							parmj = getPARM("uff", atomtypes[j])
							Xij = (parmi[0]*parmj[0])**0.5
							Dij = (parmi[1]*parmj[1])**0.5
							Aij = Dij * Xij ** 12
							Bij = -2.0 * Dij * Xij ** 6

						# Calculate scaled Electrostatic interaction (factor of 0.5 is because each interaction is counted twice)
						rij_6 = dist ** 6
						rij_12 = rij_6 ** 2
						Evw = 0.5 * VDWfactor * ((Aij / rij_12) - (Bij/rij_6))
						#print (i+1), (j+1), Evw
						VdW_term[i] = VdW_term[i] + Evw
				
				# Summary for each atom
				print "   ATOM:",str(i+1).rjust(6), atomtypes[i].rjust(10), #("%.4f" % charges[i]).rjust(10), ("%.2f" % Electrostatic_term[i]).rjust(10),
				print ("%.2f" % VdW_term[i]).rjust(10)
			print dashedline
				
				#print "\no  Total Energy:", resnumbers[residues[i][0]], resnames[residues[i][0]], "   Ees = %.5f" % Electrostatic_term[i], "   Evdw = %.5f" % VdW_term[i], "\n"
		
		## Summation over the macromolecule
		tot_Ees = sum(Electrostatic_term)
		tot_VdW = sum(VdW_term)
		print "   System Total: ","   Ees = %.5f" % tot_Ees, "   Evdw = %.5f" % tot_VdW, "\n"

if __name__ == "__main__":
	
	# Takes arguments: (1) G09 *comfile, (2) G09 *logfile, (3) *pdbfile
	# The *comfile must specify atomic connectivity - this is used to scale interactions at close range

	comfile = "none"; logfile = "none"; pdbfile = "none"
	for arg in sys.argv:
		if len(arg.split(".")) > 1:
			if arg.split(".")[1] == "com": comfile = arg.split(".")[0]

	EDA_output = calcEDA(comfile)

