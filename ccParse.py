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
#                         ccParse.py                          #
#                                                             #
#                Reads compchem job file(s)                   #
###############################################################

#Python Libraries 
import subprocess, sys, os

	
#Chemistry Libaries
from ChemUtils import *


## Check for integer when parsing ##
def is_number(s):
    try: int(s); return True
    except ValueError: return False


#Read molecule data from an input file
class getinData: 
	
	def __init__(self, file):
		if not os.path.exists(file+".com"):
			print ("\nFATAL ERROR: Input file [ %s ] does not exist"%file)

		def getJOBTYPE(self, inlines):
			for i in range(0,len(inlines)):
				if inlines[i].find("#") > -1:
					self.JOBTYPE = inlines[i].split()
			
		def getCHARGE(self, inlines):
			for i in range(0,len(inlines)):
				if inlines[i].find("#") > -1:
					#print inlines[i], inlines[i+1]
					if len(inlines[i+1].split()) == 0:
						self.CHARGE = inlines[i+4].split()[0]
						self.MULT = inlines[i+4].split()[1]
					#if len(inlines[i+2].split()) == 0:
					#	self.CHARGE = inlines[i+5].split()[0]
					#	self.MULT = inlines[i+5].split()[1]
			
		def getMEMREQ(self, inlines):
			for i in range(0,len(inlines)):
				if inlines[i].find("%mem") > -1:
					self.MEMREQ = inlines[i].split()
	
		
		def getNPROC(self, inlines):
			for i in range(0,len(inlines)):
				if inlines[i].find("%nproc") > -1:
					self.MEMREQ = inlines[i].split()
							
		
		def getATOMTYPES(self, inlines):
			self.ATOMTYPES = []
			self.LEVELTYPES = []
			for i in range(0,len(inlines)):
				if inlines[i].find("#") > -1:
					if len(inlines[i+1].split()) == 0: start = i+5
					if len(inlines[i+2].split()) == 0: start = i+6
					break
			for i in range(start,len(inlines)):
				if len(inlines[i].split()) ==0:
					break
				else:
					self.ATOMTYPES.append(inlines[i].split()[0].split("-")[0])
	
					for oniomlevel in ["H", "M", "L"]:
						if inlines[i].rfind(oniomlevel)>1:
							self.LEVELTYPES.append(inlines[i][inlines[i].rfind("H"):])
							break
		
		
		def getCONNECTIVITY(self, inlines, natoms):
			self.CONNECTIVITY = []
			self.OPTIONAL = []
			for i in range(0,len(inlines)):
				if inlines[i].find("#") > -1:
					start = i+natoms+6
					break
			
			if start < len(inlines):
				j = 1
				for i in range(start,len(inlines)):
					num = 0
					if len(inlines[i].split()) != 0:
						try:
							num = int(inlines[i].split()[0])
						except ValueError:
							num = 0
					if num == j:
						bond=[]
						neighbors=(len(inlines[i].split())-1)/2
						if neighbors!=0:
							for k in range(0,neighbors):
								bond.append((inlines[i].split()[1+2*k])+"__"+(inlines[i].split()[2+2*k]))
						self.CONNECTIVITY.append(bond)
						j = j+1
				
				if len(self.CONNECTIVITY) == natoms:
					for i in range(0, natoms): 
						for partner in self.CONNECTIVITY[i]:
							info = partner.split("__")
							nextatom = int(info[0])-1
							bondorder = float(info[1])
							nope=0
							for otherpartner in self.CONNECTIVITY[nextatom]:
								otherinfo = otherpartner.split("__")
								othernextatom = int(otherinfo[0])-1
								if othernextatom==i:
									nope=nope+1
							if nope==0:
								self.CONNECTIVITY[nextatom].append(str(i+1)+"__"+info[1])
				
				
				for i in range(start+j,len(inlines)):		
					if len(inlines[i].split()) != 0:
						self.OPTIONAL.append(inlines[i])
					
					
		def getCARTESIANS(self, inlines, natoms):
			self.CARTESIANS = []
			for i in range(0,len(inlines)):
				if inlines[i].find("#") > -1:
					start = i+5
					break
			
			for i in range(start,len(inlines)):
				if len(inlines[i].split()) == 0:
					break
				elif len(inlines[i].split()) == 4:
					self.CARTESIANS.append([float(inlines[i].split()[1]), float(inlines[i].split()[2]), float(inlines[i].split()[3])])

		def getBONDINDEX(self,inlines,natoms):
			conn=[]
			connectivity = 0
			
			for j in range(0,len(inlines)):
				if "1 2 " in inlines[j]:
					startconn = j
					connectivity  = 1
					break

			if connectivity == 1:
				for j in range(startconn,len(inlines)):
						conn.append(inlines[j])
			
			self.BONDINDEX=[]
			
			for j in range(0,natoms):
				self.BONDINDEX.append([0])
				for k in range(0,natoms):
					self.BONDINDEX[j].append(0)
			
			for j in range(0,natoms):
				if connectivity == 1:
					for bonded in conn[j].split():
						if is_number(bonded) ==True:
							if int(bonded)-1!=j:
								self.BONDINDEX[j][int(bonded)-1]=1
								self.BONDINDEX[int(bonded)-1][j]=1
		
					
		def getCONSTRAINED(self, optional):
			self.CONSTRAINED = []
			for line in optional:
				if line.find("B") > -1 and line.find("F") > -1: 
					self.CONSTRAINED.append([int(line.split(" ")[1])-1,int(line.split(" ")[2])-1])

		infile = open(file+".com","r") 
		inlines = infile.readlines()
		
		getJOBTYPE(self, inlines)
		getCHARGE(self, inlines)
		getMEMREQ(self, inlines)
		getNPROC(self, inlines)
		getATOMTYPES(self, inlines)
		self.NATOMS=len(self.ATOMTYPES)
		getCONNECTIVITY(self, inlines, self.NATOMS)
		getBONDINDEX(self,inlines,self.NATOMS)
		getCARTESIANS(self, inlines, self.NATOMS)
		getCONSTRAINED(self, self.OPTIONAL)
		
		#print "\nSuccessfully read geometry input", file
	

# Parameters specified in the params file to be used in the conformational search
class getParams:
	
	#Default search parameters
	def __init__(self, instruct):
		
		self.CSEARCH = "MCMM"
		self.MCSS="Uniform Usage Directed"
		self.EWIN=20.0
		self.COMP=10.0
		self.DEMX=41.84
		self.RJCT=0.50
		self.FMTYPE=1
		self.STEP=1000
		self.POOL=1
		self.FIXT=[]
		self.EQUI=[]
		self.RUN="/opt/mopac/MOPAC2009.exe"
		self.PROG="MOPAC2009"
		self.HSWAP=0
		self.WAIT=60
	
		#These are replaced by those in instruct file 
		try:
			n=0
			if instruct!="default":
				if not os.path.exists(instruct):
					print ("\nFATAL ERROR: Full Monte parameter file [ %s ] does not exist"%instruct)
				
				file = open(instruct,"r") 
				for line in file.readlines():
					if line.find("STEP") > -1: self.STEP=int(line.split("=")[1])
					if line.find("POOL") > -1: self.POOL=int(line.split("=")[1])
					if line.find("RJCT") > -1: self.RJCT=float(line.split("=")[1])
					if line.find("COMP") > -1: self.COMP=float(line.split("=")[1])
					if line.find("EWIN") > -1: self.EWIN=float(line.split("=")[1])
					if line.find("DEMX") > -1: self.DEMX=float(line.split("=")[1])	
					if 	line.find("FIXT") > -1: self.FIXT.append((line.split("=")[1]).rstrip('\n').lstrip())
					if 	line.find("EQUI") > -1: self.EQUI.append((line.split("=")[1]).rstrip('\n').lstrip())
					if 	line.find("HSWAP") > -1: self.HSWAP=int((line.split("=")[1]).rstrip('\n').lstrip())	
		except:
			pass	
	
		#Check that specified equivalency makes sense 
		for equi in self.EQUI:
			atomstring=[]
			equivatoms=equi.split(" ")
			for atom in equivatoms:
				if int(atom)>natoms:
					log.Fatal("\nFATAL ERROR: Full Monte parameter file [ %s ] equivalent atom numbers"%instruct)
				atomstring.append(atomtypes[int(atom)-1])
			
			for i in range(1,len(atomstring)):
				if atomstring[i]!=atomstring[i-1]:
					log.Fatal("\nFATAL ERROR: Full Monte parameter file [ %s ] equivalent atoms are different elements"%instruct)
				
		#Check that specified fixed atoms make sense
		self.FIXEDATOMS=[]
		for fix in self.FIXT:
			bondexists=0
			self.FIXEDATOMS.append([int(fix.split(" ")[0]), int(fix.split(" ")[1])])
			for bond in bondmatrix[int(fix.split(" ")[0])-1]:
				partner  = int(bond.split("__")[0])
				if partner==int(fix.split(" ")[1]):
					bondexists=bondexists+1
			if bondexists==0:
				log.Fatal("\nFATAL ERROR: Full Monte parameter file [ %s ] fixed torsion not bonded"%instruct)
			

		#print "\nSuccessfully read search parameters", instruct
	

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
                                        if outlines[i].find(" Sx") > -1:
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
			if format == "Gaussian":
				anharmonic_geom=0
				for i in range(0,len(outlines)):
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
						if anharmonic_geom==0:self.CARTESIANS.append([float(outlines[i].split()[3]),float(outlines[i].split()[4]),float(outlines[i].split()[5])])
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
					if outlines[i].find("Mulliken atomic charges:") > -1:
						self.MULLIKEN = []
						for j in range(i+2,i+natoms+2):
							self.MULLIKEN.append(float(outlines[j].split()[2]))
		
		
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
	                        	if outlines[i].find("Sx") > -1: 
						self.S2 = (float(outlines[i].split()[7]))
	
		def getENERGY(self, outlines, format):
			if format == "Gaussian":
				uff = 0
				am1 = 0
				pm3 = 0
				scf = 0
				oniom = 0
				for i in range(0,len(outlines)):
					if outlines[i].find(" UFF") > -1: uff = i
					if outlines[i] .find("AM1") > -1: am1 = i
					if outlines[i].find("PM3") > -1: pm3 = i	
					if outlines[i].find("ONIOM") > -1: oniom = i
					if outlines[i].find("SCF Done") > -1: scf = i
					if outlines[i].find("RB3LYP") > -1: self.FUNCTIONAL = "B3LYP"
				
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
		#getSOLVENERGY(self, outlines, self.FORMAT)
		getFREQS(self, outlines, self.FORMAT)
		getCPU(self, outlines, self.FORMAT)
		getATOMTYPES(self, outlines, self.FORMAT)
		if hasattr(self, "NATOMS"):
			getMULLIKEN(self, outlines, self.NATOMS, self.FORMAT)
		#getCONSTRAINED(self, outlines, self.FORMAT)
		
		#print "\nSuccessfully read geometry output", file


# Gets important data from a completed G03 outfile
def readG03out(filename,PROG): 
	#return charge, mult, atomtypes, cartesians, energy, solv, freqs if existst - gibbs, enthalpy, zpe - cpu
	semptype = 0 # Semi-Emprical
	mmechtype = 0 # Molecular Mechanics
	oniomtype = 0 # ONIOM (hybrid method)
	pcmtype = 0 # PCM solvation
	geometry=[]
	energy=0.0
	time=[0,0,0,0]
	try:
		n=0
		file = open(filename+".out","r") 
		lines = file.readlines()
		if PROG=="MOPAC2009" or PROG=="MOPAC2007":
			startgeom=0
			endgeom=0
			for line in lines:
				n=n+1
				if line.find("TOTAL ENERGY") > -1: #Get the energy (convert from eV to Hartree)
					energy=(float(line.split()[3]))	
					energy=energy*0.036749309
				if line.find("CARTESIAN COORDINATES") > -1: #Get optimized coordinates 
					startgeom=n+3
				if line.find("ATOMIC ORBITAL ELECTRON POPULATIONS") > -1: #Get optimized coordinates 
					endgeom=n-3
				if line.find("TOTAL CPU TIME") > -1:
					time=[0.0,0.0,0.0,float(line.split()[3])]
			for j in range (startgeom,endgeom):
				coordline=[]
				coordline = lines[j].split()
				#print coordline
				geometry.append([float(coordline[2]),float(coordline[3]),float(coordline[4])])	
	except:
		pass
	time=[int(time[0]), int(time[1]), int(time[2]), int(time[3])]
	return [energy, geometry, time]


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
		print file, 
		if not hasattr(fileData,"TERMINATION"): print "Warning! Job did not finish normally!"
		if hasattr(fileData, "ENERGY"): print fileData.ENERGY,
		else: print "N/A",
		if hasattr(fileData, "ZPE"): print fileData.ZPE,
		else: print "N/A", 
		if hasattr(fileData, "ENTHALPY"): print fileData.ENTHALPY,
		else: print "N/A", 
		if hasattr(fileData, "GIBBS"): print fileData.GIBBS,
		else: print "N/A", 
		if hasattr(fileData, "SOLVENERGY"): print fileData.SOLVENERGY,
		else: print "N/A", 
		if hasattr(fileData, "S2"): print fileData.S2,
                print ""
