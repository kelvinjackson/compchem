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


###############################################################
#                           FMC_1.py                          #
#              Monte Carlo Conformational Search              #
#         Dr Robert S Paton, University of Oxford 2010        #
###############################################################
#######  Written by:  Rob Paton ###############################
#######  Last modified:  Mar 20, 2013 #########################
###############################################################

 
# Python Libraries ############################################
import glob, subprocess, sys, os, random, math, tarfile
###############################################################

# Full Monte Libaries #########################################
from FMTools import *
###############################################################

if __name__ == "__main__":
	
# An input file must be specified #############################
	instruct = "default"
	interactivemode = 1
	if len(sys.argv)>1: 
		for arg in sys.argv:
			if arg == "-setup": SETUPEXE(MOPAC_EXEC)
			if arg == "-background": interactivemode = 0
		filein = sys.argv[1].split(".")[0]
		if len(sys.argv[1].split("."))>1:
			if sys.argv[1].split(".")[1] == "com": filetype = "com"
			if sys.argv[1].split(".")[1] == "pdb": filetype = "pdb"
		if len(sys.argv)>2 and sys.argv[2] != "-background": instruct = sys.argv[2]
	else: print "\nWrong number of arguments used. Correct format: FullMonte struc.com [params] \n"; sys.exit()
###############################################################
	
# Check if MOPAC.EXEC is defined ##############################
	if not os.path.exists(MOPAC_EXEC): print "\no  Mopac executable cannot be found. Rerun Full Monte with -setup as one of the arguments\n"; sys.exit()
###############################################################

# Initialize the logfile for all text output ##################
	if os.path.exists(filein+"_fm.log") and interactivemode == 1: 
		var = raw_input("\no  Log file already exists! OK to overwrite this file ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": print "   Overwriting ..."
		else: print "\nExiting\n";  sys.exit(1)
	log = FMLog(filein,"log", "fm")
###############################################################	
	
# See if there are any remaining files from a previous run ####
	if os.path.exists(filein+"_fm.tgz") and interactivemode == 1: 
		var = raw_input("\no  Tarfile already exists! OK to overwrite these structures ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": os.remove(filein+"_fm.tgz")
		else: print "\nExiting\n";  sys.exit(1)
        if len(glob.glob(filein+"*step*"))>0:
		if interactivemode == 1:
                	print glob.glob(filein+"*step*")
			var = raw_input("\no  Some optimization ouput files already exist! OK to overwrite these structures ? (Y/N) ")
                	if var.lower() == "y" or var.lower() == "": 
				for file in glob.glob(filein+"*step*"): os.remove(file)
                	else: print "\nExiting\n";  sys.exit(1)
		else:
			for file in glob.glob(filein+"*step*"): os.remove(file)
###############################################################			

# Open the structure file #####################################
	log.Write("\no  Extracting molecule from "+filein+"."+filetype+" ...")
	MOLSPEC = getinData(filein,log)
###############################################################	
	
# Open the specified parameter file for Monte Carlo parameters 
# (default values will be used if not supplied) ###############
	if instruct!="default": log.Write("\no  Extracting conformational search parameters from "+instruct+" ...")
	else: log.Write("\no  No FullMonte parameters specified! Using default values ...")
	SEARCHPARAMS = getParams(MOLSPEC, instruct,log)
###############################################################	
		
# Model Chemistry to be used ##################################
	for level in ["AM1", "PM3", "PM6", "PM7", "PM6-DH2"]:
		if SEARCHPARAMS.LEVL.upper() == level: JOB = JobSpec("Mopac")
	for level in ["UFF"]:
		if SEARCHPARAMS.LEVL.upper() == level: JOB = JobSpec("Gaussian")
	if JOB.PROGRAM != "Mopac" and JOB.PROGRAM != "Gaussian": log.Fatal("\no  "+SEARCHPARAMS.LEVL+" Level of Theory Not Yet Supported ... ")
	JOB.JOBTYPE = SEARCHPARAMS.LEVL
	log.Write("\no  Using "+JOB.JOBTYPE+" level of theory ... ")
###############################################################	
	
# For PM6, a dispersion/H-bond correction is available ########
	if JOB.PROGRAM == "Mopac" and JOB.JOBTYPE.upper() == "PM6" and interactivemode == 1:
		var = raw_input("\no  Use dispersion and H-bonding correction ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": JOB.JOBTYPE = "PM6-DH2"; log.Write("   Dispersion and H-bond correction on ... ") 
###############################################################	
	
# Solvation with CPCM #########################################
	if JOB.PROGRAM == "Gaussian" and interactivemode == 1:
		var = raw_input("\no  Use CPCM solvation ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "":
			EPS = raw_input("\n   Enter solvent name (default=diethylether) ")
			if EPS == "": EPS = "(cpcm,solvent=diethylether)"
			JOB.JOBTYPE = JOB.JOBTYPE+" scrf"+EPS; log.Write("   CPCM solvation correction on ... ")
###############################################################

# Solvation with COSMO ########################################
	if JOB.PROGRAM == "Mopac" and interactivemode == 1:
		var = raw_input("\no  Use COSMO solvation ? (Y/N) ")
		if var.lower() == "y" or var.lower() == "": 
			EPS = raw_input("\n   Enter solvent dielectric constant (default = 78.4) ")
			if EPS == "": EPS = "78.4"
			JOB.JOBTYPE = JOB.JOBTYPE+" EPS="+EPS; log.Write("   COSMO solvation correction on ... ") 
###############################################################

# MM correction for N planarity? ##############################
	if JOB.PROGRAM == "Mopac":
		for atom in MOLSPEC.ATOMTYPES:
			if atom == "N": 
				if interactivemode == 1:
					var = raw_input("\no  Use molecular mechanics correction for planar nitrogen atoms ? (Y/N) ")
					if var.lower() == "y" or var.lower() == "": JOB.JOBTYPE = JOB.JOBTYPE+" mmok"; log.Write("   MM correction on ... ") 
					else: log.Write("   No MM correction ... ") 
					break
				else:
					JOB.JOBTYPE = JOB.JOBTYPE+" mmok"; log.Write("   MM correction on ... ")
###############################################################				
	
# Check for any constraints specified #########################
	if hasattr(MOLSPEC, "CONSTRAINED"):
		JOB.CONSTRAINED = MOLSPEC.CONSTRAINED
		#print MOLSPEC.CONSTRAINED
		for const in MOLSPEC.CONSTRAINED: 
			if len(const) == 1: log.Write("\no  The Cartesian position of "+str(const[0]+1)+" will be constrained ...")
			if len(const) == 2: log.Write("\no  The distance "+str(const[0]+1)+"-"+str(const[1]+1)+" will be constrained ...")
			if len(const) == 3: log.Write("\no  The angle "+str(const[0]+1)+"-"+str(const[1]+1)+"-"+str(const[2]+1)+" will be constrained ...")
			if len(const) == 4: log.Write("\no  The dihedral "+str(const[0]+1)+"-"+str(const[1]+1)+"-"+str(const[2]+1)+"-"+str(const[3]+1)+" will be constrained ...")
		if JOB.PROGRAM == "Gaussian":
			if len(MOLSPEC.CONSTRAINED)!=0: JOB.JOBTYPE = "opt(small,modredundant,loose) "+JOB.JOBTYPE
			else: JOB.JOBTYPE = "opt(small,loose) "+JOB.JOBTYPE
###############################################################

if JOB.PROGRAM == "Gaussian": 
	JOB.JOBTYPE = "# geom=connectivity "+JOB.JOBTYPE
	JOB.NPROC = SEARCHPARAMS.POOL

# Monte Carlo or Systematic (for comparison) ##################
if interactivemode == 1:
	var = raw_input("\no  MCMM (Y) or SUMM (N) ? (Y/N) ")
	if var.lower() == "y" or var.lower() == "": SEARCHPARAMS.CSEARCH = "MCMM"
	if var.lower() == "n": SEARCHPARAMS.CSEARCH = "SUMM" 
###############################################################

# Perform an optimization of the starting geometry ############
MOLSPEC.NAME = MOLSPEC.NAME+"_step_0"
writeInput(JOB, MOLSPEC)
submitJob(JOB, MOLSPEC,log)
while isJobFinished(JOB, MOLSPEC) == 0: time.sleep(0.1)
if isJobFinished(JOB, MOLSPEC) == -1: log.Fatal("\nFATAL ERROR: Optimization of [ %s ] stuck"%file)
if isJobFinished(JOB, MOLSPEC) == 2: log.Fatal("\nFATAL ERROR: Optimization of [ %s ] failed"%file)
###############################################################	
		
# Read the output from the optimization then clean up #########
MOLSPEC.CARTESIANS =  getoutData(MOLSPEC).CARTESIANS
MOLSPEC.ENERGY =  getoutData(MOLSPEC).ENERGY
for suffix in [".com", ".mop", ".arc", ".joblog", ".errlog", ".chk"]:
	if os.path.exists(MOLSPEC.NAME+suffix): os.remove(MOLSPEC.NAME+suffix)
###############################################################
	
# Assign variable torsions, number of separate molecules and ##
# (eventually when I get round to it) rings ###################
# If number of steps is not assigned use 3^rotatable torsions #
FMVAR = Assign_Variables(MOLSPEC, SEARCHPARAMS, log)
if SEARCHPARAMS.CSEARCH == "MCMM" and SEARCHPARAMS.STEP == 0: SEARCHPARAMS.STEP = int(math.pow(3,FMVAR.MCNV))
if SEARCHPARAMS.CSEARCH == "SUMM":
	if interactivemode == 1:
		SEARCHPARAMS.ITVL = raw_input("\no  Required Interval (degrees) for systematic rotations ? ")
		if SEARCHPARAMS.ITVL == "": SEARCHPARAMS.ITVL = "60"
		SEARCHPARAMS.ITVL = int(SEARCHPARAMS.ITVL)
	interval = SEARCHPARAMS.ITVL; ninterval = 360.0/interval
	SEARCHPARAMS.STEP = int(math.pow(ninterval,FMVAR.MCNV))-1
###############################################################
		
		
# MONTE CARLO SEARCH ##########################################
start = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
asciiArt(start); Writeintro(MOLSPEC, SEARCHPARAMS, FMVAR, start, log) 

CONFSPEC = MOLSPEC
CSEARCH.NAME.append(MOLSPEC.NAME)
CSEARCH.CARTESIANS.append(MOLSPEC.CARTESIANS)
CSEARCH.TORVAL = [getTorsion(MOLSPEC)]
CSEARCH.CONNECTIVITY.append(MOLSPEC.CONNECTIVITY)
CSEARCH.ENERGY = [MOLSPEC.ENERGY]
CSEARCH.GLOBMIN = MOLSPEC.ENERGY
CSEARCH.CPU = [getoutData(MOLSPEC).CPU]
CSEARCH.ALLCPU = [getoutData(MOLSPEC).CPU]
CSEARCH.LASTFOUND = 0
CSEARCH.CLASH = [0]

print CSEARCH.STEP, SEARCHPARAMS.POOL, SEARCHPARAMS.STEP

# Stop once number of steps exceeded or no new conformers found	
while CSEARCH.STEP*SEARCHPARAMS.POOL <= SEARCHPARAMS.STEP:
	log.Write("o  STEP "+str(CSEARCH.STEP)+": Generating "+str(SEARCHPARAMS.POOL)+" structures ...")


# Setting the geometry that will be altered to generate new conformers - only relevant to MCMM
	if SEARCHPARAMS.CSEARCH == "MCMM":
		for i in range(0, CSEARCH.NSAVED):
			if CSEARCH.ENERGY[i] - CSEARCH.GLOBMIN == 0.0: startgeom = i
		
# Generate new geometries 
	for i in range(((CSEARCH.STEP-1)*SEARCHPARAMS.POOL+1),((CSEARCH.STEP)*SEARCHPARAMS.POOL)+1):
		CONFSPEC.NAME = filein+"_step_"+str(i)
		if SEARCHPARAMS.CSEARCH == "MCMM":
			if SEARCHPARAMS.MCSS == "Uniform Usage Directed":
				for j in range(0, CSEARCH.NSAVED):
					if (CSEARCH.ENERGY[j] - CSEARCH.GLOBMIN) * 2625.5 < SEARCHPARAMS.EWIN:
						if CSEARCH.USED[j] < CSEARCH.USED[startgeom]: startgeom = j
						if CSEARCH.USED[j] == CSEARCH.USED[startgeom] and CSEARCH.ENERGY[j] < CSEARCH.ENERGY[startgeom]: startgeom = j
				CSEARCH.USED[startgeom] = CSEARCH.USED[startgeom] + 1

		NBcontacts = 1

		if SEARCHPARAMS.CSEARCH == "SUMM":
			torsiontwist = [0]* len(FMVAR.TORSION)
			for j in range(0,len(FMVAR.TORSION)-1):
				#print i, math.pow(ninterval, (len(torsiontwist)-j-1))
				#print int(i)/int(math.pow(ninterval, (len(torsiontwist)-j-1))), int(i)%int(math.pow(ninterval, (len(torsiontwist)-j-1)))
				torsiontwist[j] = int(i)/int(math.pow(ninterval, (len(torsiontwist)-j-1))) * interval
				while torsiontwist[j] >= 360.0: torsiontwist[j] = torsiontwist[j] - 360.0
			torsiontwist[len(FMVAR.TORSION)-1] = int(i)%int(math.pow(ninterval, (len(torsiontwist)-j-1))) * interval
			print "   Dihedral twists in degrees:", torsiontwist

			if FMVAR.MCNV != 0: FMVAR.ADJUST = []
			for j in range(0,len(FMVAR.TORSION)): FMVAR.ADJUST.append([int(FMVAR.TORSION[j][0])+1, int(FMVAR.TORSION[j][1])+1, int(torsiontwist[j])])
			if hasattr(FMVAR, "ADJUST"):
				#print FMVAR.ADJUST
				for torsion in FMVAR.ADJUST:
					if torsion[2] != 0: CONFSPEC.CARTESIANS = AtomRot(MOLSPEC, torsion, MOLSPEC.CARTESIANS)
			NBcontacts = checkDists(CONFSPEC, SEARCHPARAMS)
			CSEARCH.CLASH.append(NBcontacts)
		
		if SEARCHPARAMS.CSEARCH == "MCMM":
			while NBcontacts > 0:
				
# The coordinates of the lowest energy, least used structure will be altered
				CONFSPEC.CARTESIANS = CSEARCH.CARTESIANS[startgeom]
				CONFSPEC.CONNECTIVITY = CSEARCH.CONNECTIVITY[startgeom]
				CONFSPEC.MMTYPES = MOLSPEC.MMTYPES
				nrandom = random.randint(FMVAR.MCNVmin, FMVAR.MCNVmax)
					
				if FMVAR.MCNV != 0:
					FMVAR.ADJUST = []
					for dihedral in random.sample(FMVAR.TORSION, nrandom): 
						FMVAR.ADJUST.append([int(dihedral[0])+1, int(dihedral[1])+1, random.randint(0,360)])
					
					if len(FMVAR.ETOZ) > 0:
						ezisomerize = random.choice([0,1])
						for dihedral in random.sample(FMVAR.ETOZ,ezisomerize): 
							#print "ETOZ",dihedral,ezisomerize
							FMVAR.ADJUST.append([int(dihedral[0]), int(dihedral[1]), 180])
				
				#Torsions in Rings - Currently Not Supported
				#if MCRI!=0:
				#	for ringtorsion in random.sample(FMVAR.RING, 1): FMVAR.ADJUST.append([int(ringtorsion[0])+1, int(ringtorsion[1])+1, random.randint(90,180)])
				#print "   Applying", len(FMVAR.ADJUST), "torsional rotations to", CSEARCH.NAME[startgeom],"..." 
				#for adj in FMVAR.ADJUST: print "   About",adj[0],"-",adj[1],"by",adj[2],"deg"
				
# Take input geometry and apply specified torsional changes
				if hasattr(FMVAR, "ADJUST"):
					print FMVAR.ADJUST
					for torsion in FMVAR.ADJUST: CONFSPEC.CARTESIANS = AtomRot(MOLSPEC, torsion, CONFSPEC.CARTESIANS)
				
# For separate molecules, alter the distances and orientations between a random number of them
				if FMVAR.NMOLS > 1:
					CONFSPEC.CARTESIANS = translateMol(FMVAR, CONFSPEC)
					CONFSPEC.CARTESIANS = rotateMol(FMVAR, CONFSPEC)

# Check for any VDW contacts smaller than specified limits
				NBcontacts = checkDists(CONFSPEC, SEARCHPARAMS)
			CSEARCH.CLASH.append(NBcontacts)

# Write input and optimize geometry
		#print CONFSPEC.CARTESIANS
		if NBcontacts == 0: writeInput(JOB, CONFSPEC); submitJob(JOB, CONFSPEC, log)

# Filter after optimization
	for i in range(((CSEARCH.STEP-1)*SEARCHPARAMS.POOL+1),((CSEARCH.STEP)*SEARCHPARAMS.POOL)+1):
		CONFSPEC.NAME = filein+"_step_"+str(i)
		
# Make sure optimization is complete
		#print i, CSEARCH.CLASH
		if CSEARCH.CLASH[i] == 0:
			while isJobFinished(JOB, CONFSPEC) == 0: time.sleep(0.1)
			
# Successful termination - extract details
			if isJobFinished(JOB, CONFSPEC) == 1:
				CONFSPEC =  getoutData(CONFSPEC)
				CSEARCH.ALLCPU.append(CONFSPEC.CPU)
				
#Check whether the molecule has isomerized - usually this is undesirable so filter out structural isomers
				concheck = checkconn(CONFSPEC, MOLSPEC, CSEARCH, SEARCHPARAMS)
				if concheck[0] == 0: CONFSPEC.CONNECTIVITY = MOLSPEC.CONNECTIVITY; isomerize = 0
				else: isomerize = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected: "+concheck[1]+concheck[2]+" has broken from "+concheck[3]+concheck[4]).ljust(50))
				
#Check whether the molecule has high energy				
				if ((CONFSPEC.ENERGY-CSEARCH.GLOBMIN)*2625.5) < SEARCHPARAMS.DEMX: toohigh = 0
				else: toohigh = 1; log.Write("   "+(CONFSPEC.NAME+" is rejected due to high energy ... ").ljust(50))
				samecheck=0

# Save or discard the optimized structure - reject if higher than global minimum by DEMX kJ/mol
				if toohigh == 0 and isomerize == 0: 
					for j in range(0, CSEARCH.NSAVED): 
						if (CONFSPEC.ENERGY-CSEARCH.GLOBMIN)*2625.5 < -0.1: break
						if abs((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5) < 0.5:
							#print str(CONFSPEC.ENERGY)+"   "+CONFSPEC.NAME+" cf "+CSEARCH.NAME[j]+" = "+str((CONFSPEC.ENERGY-CSEARCH.ENERGY[j])*2625.5)
							if checkSame(CONFSPEC, CSEARCH, SEARCHPARAMS, j) > 0 or checkSame(makemirror(CONFSPEC), CSEARCH, SEARCHPARAMS, j) > 0:
								log.Write("   "+(CONFSPEC.NAME+" is a duplicate of conformer "+CSEARCH.NAME[j]+" ... ").ljust(50))
								CSEARCH.TIMESFOUND[j] = CSEARCH.TIMESFOUND[j] + 1
								CSEARCH.NREJECT = CSEARCH.NREJECT + 1
								samecheck = samecheck + 1
								break
					
# Unique conformation with low energy! ########################			
					if samecheck == 0:
						if CONFSPEC.ENERGY < CSEARCH.GLOBMIN:
							CSEARCH.GLOBMIN = CONFSPEC.ENERGY
							log.Write("   "+(CONFSPEC.NAME+" is a new Global Minimum!").ljust(80)+("E = "+str(CSEARCH.GLOBMIN)).rjust(rightcol))
						else : log.Write("   "+(CONFSPEC.NAME+" is saved").ljust(80)+("E = "+str(CONFSPEC.ENERGY)).rjust(rightcol))
						AddConformer(CSEARCH, CONFSPEC)
						if (CONFSPEC.ENERGY-CSEARCH.GLOBMIN)*2625.5 < SEARCHPARAMS.EWIN: CSEARCH.LASTFOUND = CSEARCH.STEP*SEARCHPARAMS.POOL
###############################################################
			
# Rejection - discard #########################################
				else: CSEARCH.NREJECT = CSEARCH.NREJECT + 1
			else: 
				log.Write("\n   Unsuccessful optimization of "+CONFSPEC.NAME+".out ...")
				CSEARCH.NFAILED = CSEARCH.NFAILED + 1

			CleanAfterJob(JOB, CONFSPEC, samecheck, toohigh, isomerize)
			OrderConfs(CSEARCH, SEARCHPARAMS, start, log)
###############################################################

		else: log.Write("\n  "+CONFSPEC.NAME+" not run due to extreme steric crowding ...")
			
# End of step - update the search statistics ##################
	if (CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.NFAILED) != 0: CSEARCH.AERATE = float(CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.NREJECT-CSEARCH.NFAILED)/float(CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.NFAILED)*100.0
	else: CSEARCH.AERATE = 0.0
	if len(CSEARCH.TIMESFOUND) > 0:
		for dup in CSEARCH.TIMESFOUND:
			if dup < CSEARCH.DMIN: CSEARCH.DMIN = dup
	else: CSEARCH.DMIN = 0
###############################################################

#Tidy up the geometries and output the results ################
	if CSEARCH.STEP % 100 == 0:
		todel=[]
		#print (CSEARCH.NAME)
		for i in range(0,len(CSEARCH.NAME)):
			if ((CSEARCH.ENERGY[i] - CSEARCH.GLOBMIN)*2625.5) > SEARCHPARAMS.DEMX or (i > 199):
				#print CSEARCH.NAME[i], "earmarked for deletion"
				todel.append(i)
		if len(todel) !=0: RemoveConformer(CSEARCH, todel)
		
		WriteSummary(CSEARCH, SEARCHPARAMS, start, log)	
		CleanUp(CSEARCH, SEARCHPARAMS, filein, log)
		makeGVformat(filein, MOLSPEC, CSEARCH, SEARCHPARAMS, "fm"); makePDBformat(filein, MOLSPEC, CSEARCH, "fm")	
		
	#End of step - update step number
	CSEARCH.STEP = CSEARCH.STEP + 1
###############################################################
		
#Summary of completed Full Monte search #######################
CSEARCH.COMPLETE = 1
if SEARCHPARAMS.CSEARCH == "MCMM":
		if (CSEARCH.STEP*SEARCHPARAMS.POOL-CSEARCH.LASTFOUND) >100: log.Write("\no  Full Monte stopped finding new conformers ...")
		CSEARCH.STEP = SEARCHPARAMS.STEP
else:log.Write("\no  Full Monte completed all "+str(SEARCHPARAMS.STEP)+" steps...")
WriteSummary(CSEARCH, SEARCHPARAMS, start, log)	
CleanUp(CSEARCH, SEARCHPARAMS, filein, log)
makeGVformat(filein, MOLSPEC, CSEARCH, SEARCHPARAMS, "fm"); makePDBformat(filein, MOLSPEC, CSEARCH, "fm")
###############################################################

# Multiple Minimization with higher convergence criterion ####
multmin = 1	
if multmin == 1:
	log.Write("\no  Reoptimizing conformers with strict convergence crtieria ...")
	if JOB.PROGRAM == "Mopac": JOB.JOBTYPE = JOB.JOBTYPE+" gnorm=0.0 "
	if JOB.PROGRAM == "Gaussian": JOB.JOBTYPE = JOB.JOBTYPE.replace("loose", "")
	MultMin(CSEARCH, SEARCHPARAMS,CONFSPEC, MOLSPEC, JOB, start, log)
	#OrderConfs(CSEARCH, SEARCHPARAMS, start, log)	
	CSEARCH.GLOBMIN = CSEARCH.ENERGY[0]
###############################################################

# Final Summary of Full Monte search ##########################
WriteSummary(CSEARCH, SEARCHPARAMS, start, log)	
if os.path.isfile(filein+"_fm.tgz") == 1: os.remove(filein+"_fm.tgz")
CleanUp(CSEARCH, SEARCHPARAMS, filein, log)
makeGVformat(filein, MOLSPEC, CSEARCH, SEARCHPARAMS, "fm"); makePDBformat(filein, MOLSPEC, CSEARCH, "fm")
end = time.strftime("%Y/%m/%d %H:%M:%S", time.localtime())
asciiArt(end); log.Write(normaltermination); log.Finalize()	
###############################################################					
