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

# For reading Gaussian formatted input/output files
from ccParse import *

#Python libararies
import random, sys, os, commands, string, math

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

## Work out the atoms in the same molecules
def getMollist(bondmatrix,startatom):
	
	# The list of atoms in a molecule
	atomlist=[]
	atomlist.append(startatom)
	molecule1=[]
	nextlot=[]
	count = 0

	while count<100:
		nextlot=[]
		for atom in atomlist:
			#print atom, bondmatrix[atom]
			
			for i in range(0,len(bondmatrix[atom])):
				if bondmatrix[atom][i] == 1:
					alreadyfound = 0
					for at in atomlist:
						if i == at: alreadyfound = 1
					if alreadyfound == 0: atomlist.append(i)

		count=count+1	
	return atomlist
	

## The periodic table...
elements=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

## AMBER parameters for VdW interactions
def getAMBERPARM(type):
	AMBER_TYPE = ['C','C*','CB','CC','CN','CR','CV','CW','CA','CM','Cs','CT','F','H','H1','H2','H3','H4','H5','HA','HC','HO','HP','HS','HW','IP','K','Li','N','N2','NA','NB','N3','O','O2','OH','OS','OW','P','Rb','S','SH','FE','ZN','OM']
	AMBER_RADIUS = [1.908,1.908,1.908,1.908,1.908,1.908,1.908,1.908,1.908,1.908,3.395,1.908,1.75,0.6,1.387,1.287,1.187,1.409,1.359,1.459,1.487,0.0,1.1,0.6,0.0,1.868,2.658,1.137,1.824,1.824,1.824,1.824,1.875,1.6612,1.6612,1.721,1.6837,1.7683,2.1,2.956,2,2,1.2,1.1,1.6612]
	AMBER_WELLDEPTH = [0.086,0.086,0.086,0.086,0.086,0.086,0.086,0.086,0.086,0.086,0.0000806,0.1094,0.061,0.0157,0.0157,0.0157,0.0157,0.015,0.015,0.015,0.0157,0.0,0.0157,0.0157,0.0,0.00277,0.000328,0.0183,0.17,0.17,0.17,0.17,0.17,0.21,0.21,0.2104,0.17,0.152,0.2,0.00017,0.25,0.25,0.05,0.0125,0.21]

	for i in range(0,len(AMBER_TYPE)):
		if AMBER_TYPE[i] == type: valid = i

	return [AMBER_RADIUS[valid], AMBER_WELLDEPTH[valid]]
			

## The computation of the D3 dispersion correction
class calcEDA:
	def __init__(self, comfile, logfile, pdbfile):
		
		## Reading residue information from the pdb file
		print "\no  ASSIGNING ATOM NUMBERS BY RESIDUE: READING PDB FILE <", pdbfile,">"
		pdbData = getpdbData(pdbfile)
		print "  ", len(pdbData.RESNAME), "atoms found in the pdb file, of which", pdbData.NWAT, "are water atoms - these will be omitted ..."
		if len(pdbData.RESNAME) != len(pdbData.RESNUM): print "o  WARNING: problem assigning residue names/numbers!"; sys.exit(0)
		
		## Reading coordinate and charge information from the Gaussian09 output file
		print "\no  ASSIGNING ATOMIC CHARGES AND POSITIONS: READING G09 LOG FILE <", logfile,">"
		outData = getoutData(logfile)
		print "  ", len(outData.MULLIKEN), "atomic charges read;", len(outData.CARTESIANS), "atomic positions read"
		if len(pdbData.RESNAME) != len(outData.MULLIKEN): print "o  WARNING: PDB and G09 file have different numbers of atoms!"; sys.exit(0)
		
		## This reads from the Gaussian09 input - required for connectivity
		print "\no  ASSIGNING AMBER ATOM TYPES, ATOMIC CONNECTIVITY AND ONIOM PARTITIONS: READING G09 COM FILE <", comfile,">"
		comData = getinData(comfile)
		print "  ", len(comData.ATOMTYPES), "atomtypes read;", len(comData.LEVELTYPES), "ONIOM level types read\n"
		if len(comData.ATOMTYPES) != len(comData.LEVELTYPES): print "o  WARNING: problem assigning atom types /ONIOM levels!"; sys.exit(0)
		if len(comData.ATOMTYPES) != len(outData.MULLIKEN): print "o  WARNING: PDB and G09 file have different numbers of atoms!"; sys.exit(0)
		
		## Arrays for atoms and Cartesian coordinates ##
		connectivity = comData.BONDINDEX
		bondedlist = comData.CONNECTIVITY
		atomtypes = comData.ATOMTYPES
		natom = len(atomtypes)
		charges = outData.MULLIKEN
		leveltypes = comData.LEVELTYPES
		cartesians = outData.CARTESIANS
		
		resnames = pdbData.RESNAME
		resnumbers = pdbData.RESNUM
		residues = []
		nres = pdbData.RESNUM[-1]
		for i in range(0,nres):
			residues.append([])
			for j in range(0,len(pdbData.RESNUM)):
				if pdbData.RESNUM[j] == i+1:
					residues[i].append(j)
		
		xco=[]; yco=[]; zco=[]
		for at in cartesians:
			xco.append(at[0])
			yco.append(at[1])
			zco.append(at[2])


		# Array for interaction energies by residue
		Electrostatic_term = []
		VdW_term = []
		Eps = 332.0522173
		
		# Print energies if they exceed a certain magnitude (attractive or repulsive)
		E_printcutoff = 2.0
		
		## Compute the pairwise interaction energies
		for i in range(0,len(residues)):
			Electrostatic_term.append(0.0)
			VdW_term.append(0.0)
			## We exclude water residues
			if resnames[residues[i][0]] != "WAT":
				print "o  RESID:", resnumbers[residues[i][0]], resnames[residues[i][0]]
				#print residues[i]
				for j in residues[i]:
					print "   Atom: ", str(j+1).rjust(6),"  Type:", atomtypes[j].rjust(3),"  Layer:", leveltypes[j].rjust(2), "  Z:", "%.3f".rjust(6) % charges[j]
					for k in range(0,natom):
						## Obviously j and k cannot be equal
						if j != k:
							## Exclude interactions with water molecules
							if resnames[k] != "WAT":
								## The interaction is only evalauted across layers (i.e j and k are not in the same layer)
								if leveltypes[j] != leveltypes[k]:

									## Calculate the distance between atom j and k
									xdist = xco[j]-xco[k]; ydist = yco[j]-yco[k]; zdist = zco[j]-zco[k]
									totdist = math.pow(xdist,2)+math.pow(ydist,2)+math.pow(zdist,2)
									totdist=math.sqrt(totdist)
					
									# Interactions are by default unscaled
									Esfactor=1.0
									VDWfactor =1.0
									
									# If atoms are 1,2- 1,3- or 1,4-separated scaling factors are used
									if connectivity[j][k]==1:
										Esfactor = 0; VDWfactor = 0 #When J and K are directly bonded
										#print (j+1),"&", (k+1), "are bonded"
									else:
										for bondedjatom in bondedlist[j]:
											for bondedkatom in bondedlist[k]:
												if int(bondedjatom.split("__")[0]) == int(bondedkatom.split("__")[0]):
													Esfactor = 0; VDWfactor = 0 # When J and K are bonded to a connecting atom
													#print (j+1),"&", (k+1), "are separated by 1 atom", int(bondedjatom.split("__")[0])
													
												for bondedlatom in bondedlist[int(bondedjatom.split("__")[0])-1]:
													if int(bondedlatom.split("__")[0]) == int(bondedkatom.split("__")[0]):
														if int(bondedlatom.split("__")[0]) != int(bondedjatom.split("__")[0]):
															Esfactor = 1/1.2; VDWfactor = 1/2.0
															#print (j+1),"&", (k+1), "are separated by 2 atoms", int(bondedjatom.split("__")[0]), "&", int(bondedlatom.split("__")[0])
															
									
									# Calculate scaled Electrostatic interaction
									Ees = Esfactor*Eps*charges[j]*charges[k]/totdist
									Electrostatic_term[i] = Electrostatic_term[i] + Ees
			
									# Calculate  AMBER 6-12 VDW parameters from atom types
									parmj = getAMBERPARM(atomtypes[j])
									parmk = getAMBERPARM(atomtypes[k])
									Ejkstar = (parmj[1]*parmk[1])**0.5
									Rjkstar = parmj[0] + parmk[0]
									rjk_6 = totdist ** 6
									rjk_12 = rjk_6 ** 2
									Rjkstar_6 = Rjkstar ** 6
									Rjkstar_12 = Rjkstar_6 ** 2
									Ajk = Ejkstar * Rjkstar_12
									Bjk = 2.0 * Ejkstar * Rjkstar_6
				
									# Calculate scaled AMBER VDW interaction
									Evw = VDWfactor * ((Ajk / rjk_12) - (Bjk/rjk_6))
									VdW_term[i] = VdW_term[i] + Evw
				
									if math.fabs(Ees) > E_printcutoff or math.fabs(Evw) > E_printcutoff:
										print "   ----> ", str(k+1).rjust(6), "  Type:", atomtypes[k].rjust(3),"  Layer:", leveltypes[k].rjust(2), "  Z:", "%.3f".rjust(6) % charges[k],"  Rij: %.3f".rjust(6) % totdist, "  Ees: %.3f".rjust(6) % Ees, "  Aij:", "%.3f".rjust(6) % Ajk, "  Bij:", "%.3f".rjust(6) % Bjk,"  Evw: %.3f".rjust(6) % Evw
										
			#equation fromwilian emails
			
				print "\no  Total Energy:", resnumbers[residues[i][0]], resnames[residues[i][0]], "   Ees = %.5f" % Electrostatic_term[i], "   Evdw = %.5f" % VdW_term[i], "\n"
		
		## Summation over the macromolecule
		tot_Ees = sum(Electrostatic_term)
		tot_VdW = sum(VdW_term)
		print "   Macrmolecular Values:","   Ees = %.5f" % tot_Ees, "   Evdw = %.5f" % tot_VdW, "\n" 	

if __name__ == "__main__":
	
	# Takes arguments: (1) G09 *comfile, (2) G09 *logfile, (3) *pdbfile
	# The *comfile must specify atomic connectivity - this is used to scale interactions at close range
	# The logfile is used to extract atomic coordinates and atomic (Mulliken) charges
	# The *pdbfile must specify the residue numbers - this is used for grouping atomic contributions into residues

	if len(sys.argv) == 4:
		comfile = sys.argv[1].split(".com")[0]
		logfile = sys.argv[2].split(".log")[0]
		pdbfile = sys.argv[3].split(".pdb")[0]
	else:
		print "\nWrong number of arguments used. Correct format: ONIOM_EDA comfile logfile pdbfile\n"
		sys.exit()

	EDA_output = calcEDA(comfile, logfile, pdbfile)

