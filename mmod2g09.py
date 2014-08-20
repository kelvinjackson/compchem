#!/usr/bin/python

                  ###                     ###          ###      
                  ###                     ###          ###      
                  ###                     ###          ###      
#####b.   ####b.  ###### .d##b.  #####b.  ###  ####b.  #####b.  
### "##b     "##b ###   d##""##b ### "##b ###     "##b ### "##b 
###  ### .d###### ###   ###  ### ###  ### ### .d###### ###  ### 
### d##P ###  ### Y##b. Y##..##P ###  ### ### ###  ### ### d##P 
#####P"  "Y######  "Y### "Y##P"  ###  ### ### "Y###### #####P"  
###                                                             
###                                                             
###

###############################################################
# mmod2g09.py                                                 #
# python script to take conformers from Macromodel            # 
# conformational search within CUTOFF kJ/mol of global        #
# minimum and convert into Gaussian input format.             #
###############################################################


#Python Libraries 
import sys, os
from decimal import Decimal

#Some useful Chemistry arrays
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
	"Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

calendar=["","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]

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



#Write a new Gaussian09 input (*com) file
class writeG09input: 
	
	def __init__(self,file,Ginput,MolSpec):
		print "    Writing", file+"_"+Ginput.Append+".com\n"
		fileout = open(file+"_"+Ginput.Append+".com", "w")
		if hasattr(Ginput, "Mem"): fileout.write("%mem="+Ginput.Mem+"\n")
		if hasattr(Ginput, "Nproc"): fileout.write("%nprocshared="+Ginput.Nproc+"\n")
		if hasattr(Ginput, "Linda"): fileout.write("%nproclinda="+Ginput.Linda+"\n")
		fileout.write("# "+Ginput.Route+"\n\n")
		fileout.write(Ginput.Title+"\n\n")
		fileout.write(str(MolSpec.CHARGE)+" "+str(MolSpec.MULT)+"\n")
		for i in range(0,MolSpec.NATOMS):
			fileout.write(MolSpec.ATOMTYPES[i])
			for j in range(0,3):
				fileout.write("  "+str(Decimal(str((MolSpec.CARTESIANS[i][j])))))
			fileout.write("\n")
		fileout.write("\n")
		if len(Ginput.Optional) > 0: 
			for option in Ginput.Optional: fileout.write(option+"\n")
			fileout.write("\n")
		if len(Ginput.Freeze) > 0: 
			for frozen in Ginput.Freeze: fileout.write(frozen+"\n")
			fileout.write("\n")
		if hasattr(Ginput, "Link1"):
			fileout.write("--Link1--\n%chk="+Ginput.Link0+"\n")
			if hasattr(Ginput, "Mem"): fileout.write("%mem="+Ginput.Mem+"\n")
			if hasattr(Ginput, "Nproc"): fileout.write("%nprocshared="+Ginput.Nproc+"\n")
			if hasattr(Ginput, "Linda"): fileout.write("%nproclinda="+Ginput.Linda+"\n")
			fileout.write("# "+Ginput.Link1+"\n\n")
			fileout.write(Ginput.Title+"\n\n")
			fileout.write(str(MolSpec.CHARGE)+" "+str(MolSpec.MULT)+"\n")
			if Ginput.Link1.find("geom=check")==-1:
				for i in range(0,MolSpec.NATOMS):
					fileout.write(MolSpec.ATOMTYPES[i])
					for j in range(0,3):
						fileout.write("  "+str(Decimal(str((MolSpec.CARTESIANS[i][j])))))
					fileout.write("\n")
				fileout.write("\n")
			else: fileout.write("\n")
			if len(Ginput.Optional) > 0:
				for option in Ginput.Optional: fileout.write(option+"\n")
				fileout.write("\n")
		
		if hasattr(Ginput, "Radii"): fileout.write("radii="+Ginput.Radii+"\n\n")


#Define job specificiation for the Gaussian Calculation(s)
class createG09input: 
	
	def __init__(self,file,job):
				
		def getCutoff(self, job):
			for keyword in job:
				if 'CUTOFF' in keyword[0].upper():
					self.Cutoff = keyword[1]
			if not hasattr(self, "Cutoff"):
				self.Cutoff = "10.0"		
		
		def getAppend(self, job):
			for keyword in job:
				if 'APPEND' in keyword[0].upper():
					self.Append = keyword[1]
			if not hasattr(self, "Append"):
				self.Append = "new"		
				
				
		def getLink0(self, file, job,append):
			for keyword in job:
				if 'LINK0' in keyword[0].upper():
					if keyword[1].find(".chk") == -1:
						self.Link0 = keyword[1]+".chk"
					else: self.Link0 = keyword[1]
			if not hasattr(self, "Link0"):
				self.Link0 = file+"_"+append+".chk"		
				
		def getLink1(self, file, job):
			for keyword in job:
				if 'LINK1' in keyword[0].upper():
					self.Link1 = keyword[1]
						
	
		def getRoute(self,job):
			for keyword in job:
				if 'ROUTE' in keyword[0].upper():
					self.Route = keyword[1]
			if not hasattr(self, "Route"):
				print "\nFATAL ERROR: no Route section specified"
				sys.exit()	
	
		def getTitle(self,file, job):
			for keyword in job:
				if 'TITLE' in keyword[0].upper():
					self.Title = keyword[1]
			if not hasattr(self, "Title"): 
				self.Title = file	
				
		
		def getMem(self, job):
			for keyword in job:
				if 'MEM' in keyword[0].upper():
					self.Mem = keyword[1]
		
		
		def getNproc(self, job):
			for keyword in job:
				if 'NPROC' in keyword[0].upper():
					self.Nproc = keyword[1]
				

		def getLinda(self, job):
			for keyword in job:
				if 'LINDA' in keyword[0].upper():
					self.Linda = keyword[1]
		
	
		def getFreeze(self,job):
			self.Freeze = []
			for keyword in job:
				if 'FREEZE' in keyword[0].upper():
					self.Freeze.append(keyword[1])
			

		def getRadii(self,job):
			for keyword in job:
				if 'RADII' in keyword[0].upper():
					self.Radii = keyword[1]

	
		def getOptional(self,job):
			self.Optional = []
			for keyword in job:
				if 'OPTIONAL' in keyword[0].upper():
					self.Optional.append(keyword[1])
		

		getCutoff(self, job)
		getAppend(self, job)
		getLink0(self, file, job, self.Append)
		getRoute(self,job)
		getMem(self,job)
		getNproc(self,job)
		getLinda(self,job)
		getTitle(self, file, job)
		getFreeze(self,job)
		getRadii(self,job)
		getOptional(self,job)
		getLink1(self, file, job)
			
if __name__ == "__main__":
	
	#job specifications
	jobtype = []
	mmodfiles = []
	
	
	#input file(s)	
	#infiles = []
	
	# Takes arguments: (1) input file(s) (*maegz) (2) new job parameters
	if len(sys.argv) > 1: 
		for i in range(1,len(sys.argv)):
			if sys.argv[i][0:1] == "-" and sys.argv[i][0:3] != "--L":
				if any(sys.argv[i+1]):
					jobtype.append([sys.argv[i],sys.argv[i+1]])
			else:
				if any(sys.argv[i-1]):
					if sys.argv[i-1][0:1] != "-":
						if len(sys.argv[1].split("."))>1:
							#print sys.argv[1]
							if sys.argv[i].split(".")[1]=="maegz":
								mmodfiles.append(sys.argv[i].split(".")[0])
							
		
	else:
		print "\nWrong number of arguments used. Correct format: python mmod2g09.py maestrofile -cuttoff X.X -route \"nmr b3lyp/6-31G(d)\"\n"
		sys.exit()
	#print "1"
	#print mmodfiles
	for file in mmodfiles:
		#print "2"
		Ginput = createG09input(file, jobtype)
		qfile = open(file+"_"+Ginput.Append+".q", "w")
		print ""
		print "o   Applying an energy cutoff of", Ginput.Cutoff, "kJ/mol for the extraction of Macromodel conformers\n"
		
		#Extract archived Macromodel ouput into readable output
		
		#Parse - if energy is within CUTOFF of global minimum then save to MolSpec and create G09 input
		print "o   Extracting conformers from", file+".mae\n"
		infile = open(file+".mae","r") 
		inlines = infile.readlines()
		
		line1 =0
		line2 =0
		line3 =0
		line4 =0
		line5 =0
		
		natom=0
		nconf=0
		confenergy=[]
		
		# Define molecule
		class MOLECULE: pass
		MolSpec = MOLECULE()
				
		MolSpec.CHARGE =0
		MolSpec.MULT = 1
	
		MolSpec.ATOMTYPES=[]	
		#Find relative energy
		for i in range(0,len(inlines)):
			if inlines[i].find("f_m_ct") > -1: line1=i
			if inlines[i].find("p_m_ct") > -1: line1=i
			if inlines[i].find("r_mmod_Relative_Potential_Energy") > -1: 
				line2=i
			if inlines[i].find(":::") > -1: 
				line3=i

			if line1!=0 and line2>line1 and line3>line2:
				nconf = nconf+1
				confenergy.append(float(inlines[line2+line3-line1]))
				
				if confenergy[nconf-1] <= float(Ginput.Cutoff):
					print "    Conformer ",nconf," has relative energy", confenergy[nconf-1], "- Creating Gaussian Input"
					
					for j in range(line3,len(inlines)):
						if inlines[j].find("m_atom[") > -1: 
							line4=j
							natom = int(inlines[line4].translate(None,"m_atom[] {"))
						if inlines[j].find(":::") > -1: line5=j
						
						#print natom
						if line4!=0 and line5>line4:
							for k in range(line4,line5+20):
								if inlines[k].find("i_m_atomic_number") > -1:
									index = k-line4+4
									print index
							
							MolSpec.CARTESIANS=[]
							#print line4, line5, natom
							for k in range(line5+1,line5+natom+1):
								if len(inlines[k].split()) > 6:
									MolSpec.CARTESIANS.append([float(inlines[k].split()[2]),float(inlines[k].split()[3]),float(inlines[k].split()[4])])
								else:
									MolSpec.CARTESIANS.append([float(inlines[k].split()[1]),float(inlines[k].split()[2]),float(inlines[k].split()[3])])
								print inlines[k].split()
								print inlines[k].split()[20]
								#print periodictable[int(inlines[k].split()[20])]
								if len(inlines[k].split()) > 20:
									MolSpec.ATOMTYPES.append(periodictable[int(inlines[k].split()[21])])
								line4 =0
								line5 =0
							break						
					
					MolSpec.NATOMS =natom

					Gwrite = writeG09input(file+str(nconf), Ginput, MolSpec)
					qfile.write("g09 "+file+str(nconf)+"_"+Ginput.Append+"\n")
				else:
					print "    Conformer ",nconf," has relative energy", confenergy[nconf-1], "- Ignoring due to energy cuttof"
			
				line1 =0
				line2 =0
				line3 =0
				
				
			
	
