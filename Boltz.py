#!/usr/bin/python

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

# Copyright (c) 2012-2014 Robert Paton; All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#

###############################################################
# Boltzmann.py                                                #
# python script to perform a Boltzmann weighting over         #
# conformers computed chemical shifts                         #
#     RS Paton, Oxford 2014 - robert.paton@chem.ox.ac.uk      #
###############################################################

#Python Libraries 
import subprocess, sys, os, math
hartree_to_kcal = 2625.5/4.184
gas_const = 8.314
hartree_to_J = 2625500

##Scrape energy from an output file
class getoutData:
	
	def __init__(self, file):
		
		if not os.path.exists(file+".out"):
			if not os.path.exists(file+".log"):
				print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)
		
		
		def getFORMAT(self, outlines):
			for i in range(0,len(outlines)):
				if outlines[i].find("Gaussian") > -1: self.FORMAT = "Gaussian"; break
				if outlines[i].find("ORCA") > -1: self.FORMAT = "Orca"; break
		
		def getJOBTYPE(self, outlines, format):
			if format == "Gaussian":
				for i in range(0,len(outlines)):
					if outlines[i].find(" # ") > -1:
						self.JOBTYPE = outlines[i].lstrip(" #").rstrip("\n")
						break
				
		def getENERGY(self, outlines, format):
			if format == "Orca":
				for i in range(0,len(outlines)):
					if outlines[i].find("FINAL SINGLE POINT ENERGY") > -1 : # Get energy from HF or DFT calculation
						self.ENERGY = (float(outlines[i].split()[4]))

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
						if outlines[i].find("Energy= ") > -1 and outlines[i].find("Predicted")==-1 and outlines[i].find("Thermal")==-1: # Get energy from Semi-
							self.ENERGY = (float(outlines[i].split()[1]))
					if outlines[i].find("Total free energy in solution") > -1: 
						self.SOLVENERGY = (float(outlines[i+1].split()[7]))
		
							
		if os.path.exists(file+".out"):outfile = open(file+".out","r") 
		else: outfile = open(file+".log","r")
		
		outlines = outfile.readlines()
		
		getFORMAT(self, outlines)
		getJOBTYPE(self, outlines, self.FORMAT)
		getENERGY(self, outlines, self.FORMAT)


if __name__ == "__main__":
	
	#file(s)
	files = []
	
	# Takes arguments: (1) element of interest (2) input file(s)
	if len(sys.argv) > 2: 
		for i in range(2,len(sys.argv)):
			files.append(sys.argv[i].split(".")[0])
	else:
		print "\nWrong number of arguments used. Correct format: ccParse temperature file(s)\n"
		sys.exit()
	
	temp = float(sys.argv[1])
	energyarray = []
	relearray = []
	boltzfactor = []
	globalmin = 999.99
	boltzsum = 0
	globi = 0
	
	for file in files:
		fileData = getoutData(file)
		if hasattr(fileData, "ENERGY"): print fileData.ENERGY,
		if hasattr(fileData, "SOLVENERGY"): energyarray.append(fileData.SOLVENERGY)
                else: energyarray.append(fileData.ENERGY)
	
	if len(energyarray)==len(files):
		print "\no  BOLTZMANN FACTORS at "+str(temp)+"K"
		
		for i in range(0,len(energyarray)):
			if energyarray[i]<globalmin: globalmin = energyarray[i]; globi = i

		for i in range(0,len(energyarray)):
			boltzfactor.append(math.exp((globalmin-energyarray[i])*hartree_to_J/(gas_const*temp)))
			boltzsum = boltzsum + math.exp((globalmin-energyarray[i])*hartree_to_J/(gas_const*temp))
			relearray.append((energyarray[i]-globalmin)*hartree_to_kcal)
		
		nlessthan3 = 0
		nlessthan1 = 0
		nlessthan2 = 0
		for i in range(0,len(energyarray)):
			if (energyarray[i]-globalmin)*2625.5/4.168 < 1.0:
                                nlessthan1 = nlessthan1 + 1
			if (energyarray[i]-globalmin)*2625.5/4.168 < 2.0:
                                nlessthan2 = nlessthan2 + 1

			if (energyarray[i]-globalmin)*2625.5/4.168 < 3.0:
				nlessthan3 = nlessthan3 + 1			

		boltztot = 0
		for i in range(0,len(energyarray)): boltzfactor[i] = boltzfactor[i]/boltzsum

		#boltzfactor.sort(reverse=True)
		#print boltzfactor
		#for i in range(0,len(energyarray)):
		#	boltztot = boltztot + 100.0*boltzfactor[i]
			#print boltztot
		#	if boltztot > 95.0: break
		print "  ", "Conformer name".ljust(60), "   Energy (au)   Erel (kcal/mol)  Boltzmann factor"
		for i in range(0,len(files)):
			print "  ", files[i].ljust(60), "      %.3f" % energyarray[i], "      %.3f" % relearray[i], "        %.3f" % boltzfactor[i]

		print "\n   Most stable conformation (global minimum):", files[globi], " energy:", energyarray[i]
		print "   Conformers lying within 1 kcal/mol of the global minimum:", nlessthan1
		print "   Conformers lying within 2 kcal/mol of the global minimum:", nlessthan2
		print "   Conformers lying within 3 kcal/mol of the global minimum:", nlessthan3
		print ""
