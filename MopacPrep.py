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
#                       MopacPrep.py                          #
#                                                             #
#                   To create new *com file                   #
#                                                             #
#   Use the following keywords (case insensitive)             #
#   Append - text appended to filename			      #
#   Link0 - name of checkpoint file (default filename.chk)    #
#   Route - desired calculation type                          #
#   Mem - Max allowed RAM                                     #
#   Nproc - Number of CPUs                                    #
#   Route - desired calculation type                          #
#   Title - self-explanatory (default filename)               #
#   Charge - self-explanatory (defaults to 0)                 #
#   Mult - multiplicity (defaults to 1)                       #
#   Freeze - frozen coordinate                                #
#   Link1 - Route section for Link1                           #
#   Radii - radii defining solvent-accessible surface         #
#   More to come...                                           #
###############################################################


#Python Libraries 
import sys, os


#Chemistry Libraries
from ChemUtils import * 
from ccWrite import *
from ccParse import *


#Parse job specificiation
class getGinput: 
	
	def __init__(self,file,job):
				
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
	job = []
	
	#input file(s)
	outfiles = []
	infiles = []
	print sys.argv
	# Takes arguments: (1) input file(s) (*out) (2) new job parameters
	if len(sys.argv) > 1: 
		for i in range(1,len(sys.argv)):
			if sys.argv[i][0:1] == "-" and sys.argv[i][0:3] != "--L":
				if any(sys.argv[i+1]):
					job.append([sys.argv[i],sys.argv[i+1]])
			else:
				if any(sys.argv[i-1]):
					if sys.argv[i-1][0:1] != "-":
						if len(sys.argv[1].split("."))>1:
							print sys.argv[1]
							if sys.argv[i].split(".")[1]=="out":
								outfiles.append(sys.argv[i].split(".")[0])
							if sys.argv[i].split(".")[1]=="log":
                                                                outfiles.append(sys.argv[i].split(".")[0])
							if sys.argv[i].split(".")[1]=="com":
								infiles.append(sys.argv[i].split(".")[0])	
	else:
		print "\nWrong number of arguments used. Correct format: GaussianPrep file [new job parameters]\n"
		sys.exit()

	for file in outfiles:
		Ginput = getGinput(file, job)
		MolSpec = getoutData(file)
		Gwrite = writeMopacinput(file, Ginput, MolSpec)
	for file in infiles:
		Ginput = getGinput(file, job)
		MolSpec = getinData(file)
		MolSpec.MULT=1
		MolSpec.CHARGE=0
		Gwrite = writeMopacinput(file, Ginput, MolSpec)
	
