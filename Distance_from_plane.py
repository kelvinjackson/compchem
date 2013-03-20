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
#                 Distance_from_place.py                              #
#  A program to computes the perpendicular distance of a single atom  #
#  from a plane defined by three other atoms. This was used in the    #
#  study described in:                                                #
#  Hodgson, D. M.; Charlton, A.; Paton, R. S.; Thompson, A. S.        #
#  J. Org. Chem. 2013, 8, 1508–1518.                                  #
#######################################################################
#######  Written by:  Rob Paton #######################################
#######  Last modified:  Nov 20, 2012 #################################
#######################################################################

#Python Libraries
import subprocess, sys, os
	
#Chemistry Libaries
from ChemUtils import *

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
		#getSOLVENERGY(self, outlines, self.FORMAT)
		getFREQS(self, outlines, self.FORMAT)
		getCPU(self, outlines, self.FORMAT)
		getATOMTYPES(self, outlines, self.FORMAT)
		if hasattr(self, "NATOMS"):
			getMULLIKEN(self, outlines, self.NATOMS, self.FORMAT)
		#getCONSTRAINED(self, outlines, self.FORMAT)
		
		#print "\nSuccessfully read geometry output", file


if __name__ == "__main__":
        
        #0 file(s)
        files = []

        # Takes arguments: (1) input file(s)
        if len(sys.argv) > 1: 
                for i in range(1,2):
                        files.append(sys.argv[i].split(".")[0])
        else:
                print "\nWrong number of arguments used. Correct format: ccParse file(s)\n"
                sys.exit()

        for file in files:
                fileData = getoutData(file)
                Natom=int(sys.argv[2])-1
                atomA=int(sys.argv[3])-1
                atomB=int(sys.argv[4])-1
                atomC=int(sys.argv[5])-1
                
		print "N atom", fileData.CARTESIANS[Natom]
		print "atom A", fileData.CARTESIANS[atomA]	
		print "atom B", fileData.CARTESIANS[atomB]
		print "atom C", fileData.CARTESIANS[atomC]

		vecABx = fileData.CARTESIANS[atomA][0]-fileData.CARTESIANS[atomB][0] 
		vecABy = fileData.CARTESIANS[atomA][1]-fileData.CARTESIANS[atomB][1] 
		vecABz = fileData.CARTESIANS[atomA][2]-fileData.CARTESIANS[atomB][2] 
		#vecAB = [vecABx,vecABy,vecABz]

		vecCBx = fileData.CARTESIANS[atomC][0]-fileData.CARTESIANS[atomB][0] 
                vecCBy = fileData.CARTESIANS[atomC][1]-fileData.CARTESIANS[atomB][1] 
                vecCBz = fileData.CARTESIANS[atomC][2]-fileData.CARTESIANS[atomB][2] 
                #vecCB = [vecCBx,vecCBy,vecCBz]

		crossx = vecABy*vecCBz - vecABz*vecCBy
		crossy = vecABz*vecCBx - vecABx*vecCBz
		crossz = vecABx*vecCBy - vecABy*vecCBx

		crossmag = crossx**2 + crossy**2 + crossz**2
		crossmag = crossmag**0.5

		crossx = crossx/crossmag
		crossy = crossy/crossmag
		crossz = crossz/crossmag
		cross = [crossx,crossy,crossz]
		#print cross
		
		vecANx = fileData.CARTESIANS[atomA][0]-fileData.CARTESIANS[Natom][0] 
                vecANy = fileData.CARTESIANS[atomA][1]-fileData.CARTESIANS[Natom][1] 
                vecANz = fileData.CARTESIANS[atomA][2]-fileData.CARTESIANS[Natom][2] 
                #vecAN = [vecANx,vecANy,vecANz]

		dist = crossx*vecANx + crossy*vecANy + crossz*vecANz

		print dist	
