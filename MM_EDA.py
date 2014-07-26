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
#                        ONIOM_EDA.py                                 #
#     Wilian Cortopassi & Robert Paton, University of Oxford, 2014    #
#                                                                     #
#    For Gaussian formatted input/output file this program will       #
#    compute pairwise electrostatic and VdW interactions              #
#    Point charges are (at present Mulliken) taken from the G09 output#
#    1,2- and 1,3- interactions are set to zero, with 1,4-            #
#    interactions scaled according to Hirao JPCA                      #
#    As per Hirao, interactions are only evaluated between atoms in   #
#    different layers (i.e. High and Low)                             #       
#    The connectivity is required in the G09 input                    #
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
	def __init__(self, comfile, logfile, pdbfile):
		
		## Reading residue information from the pdb file
		if pdbfile != "none":
			print "\no  SCRAPING PDB FILE <", pdbfile,">"
			pdbData = getpdbData(pdbfile)
			if hasattr(pdbData, "ATOMID"): print "   Read", len(pdbData.ATOMID), "atoms",
			if hasattr(pdbData, "NMOL"): print "   Read", pdbData.NMOL, "separate molecules (including waters)",; nmols = pdbData.NMOL; mols = pdbData.MOLS
				#and", pdbData.NWAT, "are water atoms"
			if hasattr(pdbData, "RESNUM"): print "   Read", pdbData.RESNUM, "separate residues (including waters)",; nres = pdbData.RESNUM; res = pdbData.RES; resname = pdbData.RESNAME
				
		
		## Reading Gaussian09 input - connectivity information
		if comfile != "none":
			print "\no  SCRAPING G09 COM FILE <", comfile,">"
			comData = getinData(comfile)
			if hasattr(comData, "ATOMTYPES"): print "   Read", len(comData.ATOMTYPES), "atomtypes",; atomtypes = comData.ATOMTYPES
			if hasattr(comData, "LEVELTYPES"): print "   Read", len(comData.LEVELTYPES), "ONIOM level types",; qmmm_level = comData.LEVELTYPES; nlevels = comData.NLEVELS
			if hasattr(comData, "CONNECTIVITY"): print "   Read", len(comData.CONNECTIVITY), "Atomic Connectivity Information",; bondedlist = comData.CONNECTIVITY; #connectivity = comData.BONDINDEX
			#if len(comData.ATOMTYPES) != len(comData.LEVELTYPES): print "o  WARNING: problem assigning atom types /ONIOM levels!"; sys.exit(0)
			#if len(comData.ATOMTYPES) != len(outData.MULLIKEN): print "o  WARNING: PDB and G09 file have different numbers of atoms!"; sys.exit(0)
		
	
		## Reading coordinate and charge information from the Gaussian09 output file
		if logfile != "none":
			print "\no  SCRAPING G09 LOG FILE <", logfile,">"
			outData = getoutData(logfile)
			if hasattr(outData, "CARTESIANS"): print "   Read", len(outData.CARTESIANS), "Atomic Cartesian coordinates",; cartesians = outData.CARTESIANS; natom = len(cartesians)
			if hasattr(outData, "MULLIKEN"): print "  ", len(outData.MULLIKEN), "Mulliken atomic charges",; mulliken = outData.MULLIKEN
			if hasattr(outData, "APT"): print "  ", len(outData.APT), "APT atomic charges",; apt = outData.APT
	
			if pdbfile != "none" and len(pdbData.ATOMID) != natom: print "\no  WARNING: PDB and G09 file have different numbers of atoms!"; sys.exit(0)
			
			try: bondedlist
			except NameError:
				outData.CONNECTIVITY = []; #outData.BONDINDEX = []
				for i in range(0,len(outData.CARTESIANS)):
					#print i,
					outData.CONNECTIVITY.append([]); #outData.BONDINDEX.append([])
					for j in range(0,len(outData.CARTESIANS)):
						#outData.BONDINDEX[i].append(0)
						if numpy.linalg.norm(numpy.array(outData.CARTESIANS[i])-numpy.array(outData.CARTESIANS[j])) < 0.5*(bondiRadius(atomicnumber(outData.ATOMTYPES[i])) + bondiRadius(atomicnumber(outData.ATOMTYPES[j]))):
							if i != j: outData.CONNECTIVITY[i].append(str(j+1)+"__1.0"); #outData.BONDINDEX[i][j] = 1

				print "\n   Generated atomic connectivity information from Cartesian coordinates"; bondedlist = outData.CONNECTIVITY; #connectivity = outData.BONDINDEX


			try: atomtypes
			except NameError:
				atomtypes = []
				for i in range(0,len(outData.ATOMICTYPES)):
					if outData.ATOMICTYPES[i] == 0 : atomtypes.append(outData.ATOMTYPES[i])
					else: atomtypes.append(at_to_string(outData.ATOMICTYPES[i]))
		
		
		## Define number of ONIOM layers
		try: qmmm_level
		except NameError: nlevels = 1
		print "\n   There are", nlevels, "subregions (QM/MM) present:"
		if nlevels > 1:
			high_level = []; medium_level = []; low_level = []
			for i in range(0,natom):
				if qmmm_level[i] == "H": high_level.append(i)
				if qmmm_level[i] == "M": medium_level.append(i)
				if qmmm_level[i] == "L": low_level.append(i)
		
				
		## Define residues
		try: nres
		except NameError: nres = 0
		#if nres > 0:
		#	for i in range(0,nres):
		#		print i, resname[res[i][0]], res[i]
				
		
		## Define number of separate molecules
		try: nmols
		except NameError:
			molA = getMollist(connectivity,0,natom)
			if len(molA) < natom:
				mols = []
				mols.append(sorted(molA))
				for i in range(1,natom):
					found = 0
					#print "   checking atom", i
					for j in range(0,len(mols)):
						#print mols[j]
						for atom in mols[j]:
							if atom == i:
								#print "previously found", i;
								found = 1

					if found == 0:
						#print "not previously found, saving"
						nextMol = getMollist(connectivity,i,natom)
						mols.append(sorted(nextMol))

				nmols = len(mols)

		print "\n   There are", nmols, "molecules present:"
		for i in range(0,len(mols)): print "   Molecule", (i+1), "containing",len(mols[i]), "atoms"
		
		
		## Define Charges to be used by user input
		if interactivemode == 1:
			var = raw_input("\no  Define Charge Model to be Used: (1) Mulliken; (2) APT; (3) QEQ; (4) Manually define ... ")
			try: charge_model = int(var)
			except OSError, e: print >>sys.stderr; sys.exit()
				
			if charge_model == 1:
				try: charges = mulliken
				except NameError: print "problem!"
			if charge_model == 2:
				try: charges = apt
				except NameError: print "problem!"
			if charge_model == 4: print "how we gonna define?" ## to do
			if charge_model == 3: print "using chornsby's model" ## to do
	
			var = raw_input("\no  Define MM atom types to be Used: (1) AMBER; (2) UFF; (3) None ... ")
			try: mm_type = int(var)
			except OSError, e: print >>sys.stderr; sys.exit()

			var = raw_input("\no  Perform Energy Decomposition Analysis (EDA) by: (1) Atom; (2) Molecule; (3) Residue; (4) User-defined ... ")
			try: eda_model = int(var)
			except OSError, e: print >>sys.stderr; sys.exit()
		
		
		# Array for interaction energies by residue
		Electrostatic_term = []
		VdW_term = []


		## Interatomic EDA
		## Compute the pairwise interaction energies between atoms i and j
		if eda_model == 1:
			print "\n                          ", ("Charge").rjust(10), ("E(el)").rjust(10), ("E(vdw)").rjust(10)
			for i in range(0,natom):
				VdW_term.append(0.0); Electrostatic_term.append(0.0)
				for j in range(0,natom):
					if i != j:
						dist = numpy.linalg.norm(numpy.array(outData.CARTESIANS[i])-numpy.array(outData.CARTESIANS[j]))
						
						# If atoms are 1,2- 1,3- or 1,4-separated, scaling factors are used
						Esfactor=1.0; VDWfactor =1.0
						
						if connectivity[i][j]==1: Esfactor = 0; VDWfactor = 0 #When I and J are directly bonded
	
						else:
							for bondediatom in bondedlist[i]:
								for bondedjatom in bondedlist[j]:
									if int(bondediatom.split("__")[0]) == int(bondedjatom.split("__")[0]):
										Esfactor = 0; VDWfactor = 0 # When I and J are bonded to a connecting atom
									
									for bondedkatom in bondedlist[int(bondediatom.split("__")[0])-1]:
										if int(bondedkatom.split("__")[0]) == int(bondedjatom.split("__")[0]):
											if int(bondedkatom.split("__")[0]) != int(bondediatom.split("__")[0]):
												Esfactor = 1/1.2; VDWfactor = 1/2.0 # When I and J are bonded to 2 connecting atoms
						
						# Calculate scaled Electrostatic interaction (factor of 0.5 is because each interaction is counted twice)
						Ees = 0.5 * Esfactor*Eps*charges[i]*charges[j]/dist
						#print str(i+1).rjust(6), str(j+1).rjust(6), str(Esfactor).rjust(6), str(charges[i]).rjust(10), str(charges[j]).rjust(10), str(Ees).rjust(10)
						Electrostatic_term[i] = Electrostatic_term[i] + Ees

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
						VdW_term[i] = VdW_term[i] + Evw
				
				# Summary for each atom
				print "   ATOM:",str(i+1).rjust(6), atomtypes[i].rjust(10), ("%.4f" % charges[i]).rjust(10), ("%.2f" % Electrostatic_term[i]).rjust(10), ("%.2f" % VdW_term[i]).rjust(10)
			print dashedline


		## Intermolecular EDA
		## Compute the pairwise interaction energies between atoms i and j only if they are in separate molecules
		if eda_model == 2:
			for i in range(0,natom): VdW_term.append(0.0); Electrostatic_term.append(0.0)
			print "\n                          ", ("Charge").rjust(10), ("E(el)").rjust(10), ("E(vdw)").rjust(10)
			for n in range(0,len(mols)):
				Mol_Ees_term = 0.0; Mol_charge = 0.0; Mol_VdW_term = 0.0
				for i in mols[n]:
					Mol_charge = Mol_charge + charges[i]
					for o in range(0,len(mols)):
						if o != n:
							for j in mols[o]:
								
								dist = numpy.linalg.norm(numpy.array(outData.CARTESIANS[i])-numpy.array(outData.CARTESIANS[j]))
								Esfactor=1.0; VDWfactor =1.0
																
								# Calculate scaled Electrostatic interaction (factor of 0.5 is because each interaction is counted twice)
								Ees = 0.5 * Esfactor*Eps*charges[i]*charges[j]/dist
								#print str(i+1).rjust(6), str(j+1).rjust(6), str(Esfactor).rjust(6), str(charges[i]).rjust(10), str(charges[j]).rjust(10), str(Ees).rjust(10)
								Electrostatic_term[i] = Electrostatic_term[i] + Ees
								Mol_Ees_term = Mol_Ees_term + Ees
											
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
								#print i, j, Aij, Bij, rij_12, rij_6, Evw
								VdW_term[i] = VdW_term[i] + Evw
								Mol_VdW_term = Mol_VdW_term + Evw

					print "   ATOM:",str(i+1).rjust(6), atomtypes[i].rjust(10), ("%.4f" % charges[i]).rjust(10), ("%.2f" % Electrostatic_term[i]).rjust(10), ("%.2f" % VdW_term[i]).rjust(10)
				print dashedline
				print "   Mol. "+str(n+1)+" Total:".rjust(17), ("%.4f" % Mol_charge).rjust(10), ("%.2f" % Mol_Ees_term).rjust(10), ("%.2f" % Mol_VdW_term).rjust(10),
				print "\n", dashedline, "\n"
				

		## Inter-residue EDA
		## Compute the pairwise interaction energies between atoms i and j only if there are in separate residues
		if eda_model == 3:
			for i in range(0,natom): VdW_term.append(0.0); Electrostatic_term.append(0.0)
			print "\n                          ", ("Charge").rjust(10), ("E(el)").rjust(10), ("E(vdw)").rjust(10)
			for n in range(0,nres):
				Res_Ees_term = 0.0; Res_charge = 0.0; Res_VdW_term = 0.0
				for i in res[n]:
					Res_charge = Res_charge + charges[i]
						
					for o in range(0,nres):
						if o != n:
							for j in res[o]:
								if qmmm_level[i] != qmmm_level[j]:
									bonded = 0
									for bond in (bondedlist[i]):
										#print i "is bonded to", bond
										if int(bond.split("__")[0])-1 == j: bonded = 1

									dist = numpy.linalg.norm(numpy.array(outData.CARTESIANS[i])-numpy.array(outData.CARTESIANS[j]))
					
									# If atoms are 1,2- 1,3- or 1,4-separated, scaling factors are used
									Esfactor=1.0; VDWfactor =1.0
									
									if bonded==1: Esfactor = 0; VDWfactor = 0 #When I and J are directly bonded
									
									else:
										for bondediatom in bondedlist[i]:
											for bondedjatom in bondedlist[j]:
												if int(bondediatom.split("__")[0]) == int(bondedjatom.split("__")[0]):
													Esfactor = 0; VDWfactor = 0 # When I and J are bonded to a connecting atom
												
												for bondedkatom in bondedlist[int(bondediatom.split("__")[0])-1]:
													if int(bondedkatom.split("__")[0]) == int(bondedjatom.split("__")[0]):
														if int(bondedkatom.split("__")[0]) != int(bondediatom.split("__")[0]):
															Esfactor = 1/1.2; VDWfactor = 1/2.0 # When I and J are bonded to 2 connecting atoms
									
									# Calculate scaled Electrostatic interaction (factor of 0.5 is because each interaction is counted twice)
									Ees = 0.5 * Esfactor*Eps*charges[i]*charges[j]/dist
									#print str(i+1).rjust(6), str(j+1).rjust(6), str(Esfactor).rjust(6), str(charges[i]).rjust(10), str(charges[j]).rjust(10), str(Ees).rjust(10)
									Electrostatic_term[i] = Electrostatic_term[i] + Ees
									Res_Ees_term = Res_Ees_term + Ees
									
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
									VdW_term[i] = VdW_term[i] + Evw
									Res_VdW_term = Res_VdW_term + Evw
									
					print "   ATOM:",str(i+1).rjust(6), atomtypes[i].rjust(10), ("%.4f" % charges[i]).rjust(10), ("%.2f" % Electrostatic_term[i]).rjust(10), ("%.2f" % VdW_term[i]).rjust(10)
				print dashedline
				print "   Res. "+str(n+1).rjust(6)+" Total:     ", ("%.4f" % Res_charge).rjust(10), ("%.2f" % Res_Ees_term).rjust(10), ("%.2f" % Res_VdW_term).rjust(10),
				print "\n", dashedline, "\n"
			
				
				#print "\no  Total Energy:", resnumbers[residues[i][0]], resnames[residues[i][0]], "   Ees = %.5f" % Electrostatic_term[i], "   Evdw = %.5f" % VdW_term[i], "\n"
		
		## Summation over the macromolecule
		tot_Ees = sum(Electrostatic_term)
		tot_VdW = sum(VdW_term)
		print "   System Total: ","   Ees = %.5f" % tot_Ees, "   Evdw = %.5f" % tot_VdW, "\n"

if __name__ == "__main__":
	
	# Takes arguments: (1) G09 *comfile, (2) G09 *logfile, (3) *pdbfile
	# The *comfile must specify atomic connectivity - this is used to scale interactions at close range
	# The logfile is used to extract atomic coordinates and atomic (Mulliken) charges
	# The *pdbfile must specify the residue numbers - this is used for grouping atomic contributions into residues

	comfile = "none"; logfile = "none"; pdbfile = "none"
	for arg in sys.argv:
		if len(arg.split(".")) > 1:
			if arg.split(".")[1] == "com": comfile = arg.split(".")[0]
			if arg.split(".")[1] == "log" or arg.split(".")[1] == "out": logfile = arg.split(".")[0]
			if arg.split(".")[1] == "pdb": pdbfile = arg.split(".")[0]


	EDA_output = calcEDA(comfile, logfile, pdbfile)

