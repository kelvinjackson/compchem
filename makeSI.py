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

###############################################################
#                          makeSI.py                          #
#                                                             #
#        Process G09 output into Supporting Info format       #
###############################################################

#Python Libraries 
import subprocess, sys, os, commands, string, math
from decimal import Decimal

#Some useful arrays
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
	"Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

atomicmass = [1.008, 4.003, 6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07, 35.45, 39.95, 39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.84, 58.93, 58.69, 
	63.55, 65.39, 69.72, 72.61, 74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 99.0, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29]

def elementID(massno):
		if massno < len(periodictable): return periodictable[massno]	
		else: return "XX"

def atomicnumber(element):
	atomicno = 0
	for i in range(0,len(periodictable)):
		if element == periodictable[i]: atomicno = i
	return atomicno
	
## Check for integer when parsing ##
def is_number(s):
    try: int(s); return True
    except ValueError: return False


class writeXYZ:
	def __init__(self,file,MolSpec,Estring,ZPEstring,Hstring,Gstring):
		print "\nWriting", file+".xyz\n"
		fileout = open(file+".xyz", "w")
		fileout.write(str(MolSpec.NATOMS)+"\n")
		fileout.write(Estring+" "+ZPEstring+" "+Hstring+" "+Gstring+" \n")
		for i in range(0,MolSpec.NATOMS):
			fileout.write(MolSpec.ATOMTYPES[i])
			for j in range(0,3):
				fileout.write("  "+str(Decimal(str((MolSpec.CARTESIANS[i][j])))))
			fileout.write("\n")

#Read molecule data from an output file
class getoutData:
	
	def __init__(self, file):
	
		if not os.path.exists(file+".out"):
			if not os.path.exists(file+".log"):
                        	print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)

		
		def getFORMAT(self, outlines):
			for i in range(0,len(outlines)):
				if outlines[i].find("Gaussian") > -1: self.FORMAT = "Gaussian"; break
				if outlines[i].find("MOPAC") > -1: self.FORMAT = "Mopac"; break
					
					
		def getSPIN(self, outlines, format):
                        if format == "Gaussian":
                                for i in range(0,len(outlines)):
                                        if outlines[i].find("<Sx>") > -1:
                                                print outlines[i].split()
                                                self.S2 = (float(outlines[i].split()[7]))

		def getJOBTYPE(self, outlines, format):
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find(" # ") > -1:
						self.JOBTYPE = outlines[i].lstrip(" #").rstrip("\n")
						break

		def getTERMINATION(self, outlines,format):
                	if format == "Gaussian":
      				for i in range(0,len(outlines)):
                                	if outlines[i].find("Normal termination") > -1:
                                        	self.TERMINATION = "normal"

		
		def getCHARGE(self, outlines, format):
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("Charge = ") > -1:
						self.CHARGE = int(outlines[i].split()[2])
						self.MULT = int(outlines[i].split()[5].rstrip("\n"))
						break
			if format == "Mopac":
				self.CHARGE = 0
				#ideally add up all the atomic charges here?
                                self.MULT =  1 
		
								
		def getATOMTYPES(self, outlines, format):
			self.ATOMTYPES = []
			self.CARTESIANS = []
			self.ATOMICTYPES = []
			if format == "Gaussian":
				anharmonic_geom=0
				for i in range(0,len(outlines)):
					if outlines[i].find("Input orientation") > -1:
						standor = i
					if outlines[i].find("Standard orientation") > -1:
						standor = i
					if outlines[i].find("Vib.Av.Geom.") > -1:
						standor = i
						anharmonic_geom=1
					if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1:
						self.NATOMS = i-standor-6
				try: standor
				except NameError: pass
				else:
					for i in range (standor+5,standor+5+self.NATOMS):
						self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
						self.ATOMICTYPES.append(int(outlines[i].split()[2]))
						
						if anharmonic_geom==0:
							if len(outlines[i].split()) > 5: self.CARTESIANS.append([float(outlines[i].split()[3]),float(outlines[i].split()[4]),float(outlines[i].split()[5])])
							else: self.CARTESIANS.append([float(outlines[i].split()[2]),float(outlines[i].split()[3]),float(outlines[i].split()[4])])
						if anharmonic_geom==1:self.CARTESIANS.append([float(outlines[i].split()[2]),float(outlines[i].split()[3]),float(outlines[i].split()[4])])

			if format == "Mopac":
                		for i in range(0,len(outlines)):
                    			#if outlines[i].find("TOTAL ENERGY") > -1: #Get the energy (convert from eV to Hartree)
                    			#   energy=(float(line.split()[3]))	
                    			#    energy=energy*0.036749309
                    			if outlines[i].find("CARTESIAN COORDINATES") > -1: startgeom = i+4
                    			if outlines[i].find("ATOMIC ORBITAL ELECTRON POPULATIONS") > -1: endgeom = i-2; break
                    			#if outlines[i].find("TOTAL CPU TIME") > -1:
                    			#    time=[0.0,0.0,0.0,float(line.split()[3])]
                		self.NATOMS = endgeom - startgeom
				for i in range (startgeom,endgeom):
					self.ATOMTYPES.append(outlines[i].split()[1])
                    			self.CARTESIANS.append([float(outlines[i].split()[2]),float(outlines[i].split()[3]),float(outlines[i].split()[4])])	
			
														
		def getFREQS(self, outlines, format):
			self.FREQS = []
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("Frequencies") > -1:
						self.FREQS.append(float(outlines[i].split()[2]))
						if len(outlines[i].split()) > 3: self.FREQS.append(float(outlines[i].split()[3]))
						if len(outlines[i].split()) > 4: self.FREQS.append(float(outlines[i].split()[4]))
				if len(self.FREQS) > 0:
					for i in range(0,len(outlines)):
						if outlines[i].find("Zero-point correction") > -1: self.ZPE = float(outlines[i].split()[2])
						if outlines[i].find("thermal Enthalpies") > -1: self.ENTHALPY = float(outlines[i].split()[6])
						if outlines[i].find("thermal Free Energies") > -1: self.GIBBS = float(outlines[i].split()[7])


		def getMULLIKEN(self, outlines, natoms, format):
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("Mulliken charges:") > -1 or outlines[i].find("Mulliken charges and spin densities:") > -1:
						
						self.MULLIKEN = []
						for j in range(i+2,i+natoms+2):
							self.MULLIKEN.append(float(outlines[j].split()[2]))
		
		def getAPT(self, outlines, natoms, format):
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find("APT charges:") > -1:
						self.APT = []
						for j in range(i+2,i+natoms+2):
							self.APT.append(float(outlines[j].split()[2]))


		def getCPU(self, outlines, format):
			days = 0
			hours = 0
			mins = 0
			secs = 0
			if format == "Gaussian":
					for i in range(0,len(outlines)):
						if outlines[i].find("Job cpu time") > -1:
							days = days + int(outlines[i].split()[3])
							hours = hours + int(outlines[i].split()[5])
							mins = mins + int(outlines[i].split()[7])
							secs = secs + int(float(outlines[i].split()[9]))
			self.CPU=[days,hours,mins,secs]
		
		def getSPIN(self, outlines, format):
                	if format == "Gaussian":	
				for i in range(0,len(outlines)):
	                        	if outlines[i].find("<Sx>") > -1: 
						self.S2 = (float(outlines[i].split()[7]))
	
		def getENERGY(self, outlines, format):
			if format == "Gaussian":
				uff = 0
				am1 = 0
				pm3 = 0
				scf = 0
				oniom = 0
				for i in range(0,len(outlines)):
					if outlines[i].find("Standard basis:") > -1: self.BASIS = outlines[i].split()[2]
					if outlines[i].find(" UFF") > -1: uff = i
					if outlines[i] .find("AM1") > -1: am1 = i
					if outlines[i].find("PM3") > -1: pm3 = i	
					if outlines[i].find("ONIOM") > -1: oniom = i
					if outlines[i].find("SCF Done") > -1: scf = i
					if outlines[i].find("(RB3LYP)") > -1 or outlines[i].find("(UB3LYP)") > -1: self.FUNCTIONAL = "B3LYP"
					if outlines[i].find("(RB-P86)") > -1: self.FUNCTIONAL = "BP86"
					if outlines[i].find("(RB2PLYP)") > -1: self.FUNCTIONAL = "B2PLYP"
					if outlines[i].find("(RM06)") > -1: self.FUNCTIONAL = "M06"
					if outlines[i].find("(RM062X)") > -1: self.FUNCTIONAL = "M06-2X"
					if outlines[i].find("(RM06L)") > -1: self.FUNCTIONAL = "M06L"
					if outlines[i].find("(RB97D)") > -1: self.FUNCTIONAL = "B97D"
					if outlines[i].find("(RwB97XD)") > -1: self.FUNCTIONAL = "wB97XD"
				
				calctype = [uff,am1,pm3,oniom,scf]
				for i in range(0,len(outlines)):
					if scf == max(calctype) and outlines[i].find("SCF Done") > -1 and outlines[i].find("Initial convergence to 1.0D-05 achieved")==-1: # Get energy from HF or DFT calculation
						self.ENERGY = (float(outlines[i].split()[4]))
					if oniom == max(calctype) and outlines[i].find("ONIOM: extrapolated energy") > -1: # Get energy from ONIOM calculation
						self.ENERGY = (float(outlines[i].split()[4]))		
					if pm3 == max(calctype) or am1 == max(calctype) or uff == max(calctype):
						if outlines[i].find("Energy= ") > -1 and outlines[i].find("Predicted")==-1 and outlines[i].find("Thermal")==-1: # Get energy from Semi-empirical or Molecular Mechanics calculation
							self.ENERGY = (float(outlines[i].split()[1]))
					if outlines[i].find("Total free energy in solution") > -1: 
						self.SOLVENERGY = (float(outlines[i+1].split()[7]))

				
		if os.path.exists(file+".out"):outfile = open(file+".out","r") 
		else: outfile = open(file+".log","r")
		
		outlines = outfile.readlines()
		
		getFORMAT(self, outlines)
		getJOBTYPE(self, outlines, self.FORMAT)
		getTERMINATION(self, outlines,self.FORMAT)
		getCHARGE(self, outlines, self.FORMAT)
		getENERGY(self, outlines, self.FORMAT)
		getSPIN(self, outlines, self.FORMAT)
		getFREQS(self, outlines, self.FORMAT)
		getCPU(self, outlines, self.FORMAT)
		getATOMTYPES(self, outlines, self.FORMAT)
		getMULLIKEN(self, outlines, self.NATOMS, self.FORMAT)
		getAPT(self, outlines, self.NATOMS, self.FORMAT)
		#getCONSTRAINED(self, outlines, self.FORMAT)

if __name__ == "__main__":
	
	#0 file(s)
	files = []

	# Takes arguments: (1) input file(s)
	if len(sys.argv) > 1: 
		for i in range(1,len(sys.argv)):
			files.append(sys.argv[i].split(".")[0])
	else:
		print "\nWrong number of arguments used. Correct format: ccParse file(s)\n"
		sys.exit()

	for file in files:
		fileData = getoutData(file)
		natoms = len(fileData.ATOMTYPES)
		functional = fileData.FUNCTIONAL; basis = fileData.BASIS
		if not hasattr(fileData,"TERMINATION"): print "Warning! Job did not finish normally!"
		if hasattr(fileData, "ENERGY"): Estring = "E("+functional+"/"+basis+") = "+str(fileData.ENERGY)
		else: Estring = ""
		if hasattr(fileData, "ZPE"): ZPEstring = "ZPE("+functional+"/"+basis+") = "+str(fileData.ZPE)
		else: ZPEstring = ""
		if hasattr(fileData, "ENTHALPY"): Hstring = "H("+functional+"/"+basis+") = "+str(fileData.ENTHALPY)
		else: Hstring = ""
		if hasattr(fileData, "GIBBS"): Gstring = "G("+functional+"/"+basis+") = "+str(fileData.GIBBS)
		else: Gstring = ""
		IMF=[]
		if hasattr(fileData, "FREQS"): 
			for i in fileData.FREQS:
				if i <=-40: Fstring = "Imaginary frequency = "+str(i);IMF.append(i)#Set imaginary frequency cutoff
			
		if len(IMF)==0: Fstring =""
		print file+":"
		print Estring
		print ZPEstring
		print Hstring
		print Gstring
		print Fstring
		print ""
		for i in range(0,fileData.NATOMS):
			print fileData.ATOMTYPES[i],
			for j in range(0,3):
				cart = (fileData.CARTESIANS[i][j])
				print "", "%.6f" % cart,
			print ""

		writeXYZ(file, fileData, Estring, ZPEstring, Hstring, Gstring)

