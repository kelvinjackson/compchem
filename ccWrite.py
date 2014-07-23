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
#                         ccWrite.py                          #
#                                                             #
#               Writes compchem job file(s)                   #
###############################################################

from decimal import Decimal


#Create Gaussian com file
class writeGinput: 

############################################################################
#                              Gaussian Format                             #
#   %chk = $Link0                                                          #
#   (%mem = $Mem)                                                          #
#   (%nproc = $Nproc)                                                      #
#   Route - desired calculation type                                       #
#                                                                          #
#   $Charge $Mult                                                          #
#   Molecule specification                                                 #
#                                                                          #
#   ($Optional)                                                            #
#                                                                          #
#   ($Freeze)                                                              #
#                                                                          #
#   (radii = $Radii)                                                       #
#                                                                          #
############################################################################
	
	def __init__(self,file,Ginput,MolSpec):
		print "\nWriting", file+"_"+Ginput.Append+".com\n"
		fileout = open(file+"_"+Ginput.Append+".com", "w")
		fileout.write("%chk="+Ginput.Link0+"\n")
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
		
#Create Mopac mop file
class writeMopacinput:

############################################################################
#                              Mopac Format                               #
#   %chk = $Link0                                                          #
#   (%mem = $Mem)                                                          #
#   (%nproc = $Nproc)                                                      #
#   Route - desired calculation type                                       #
#                                                                          #
#   $Charge $Mult                                                          #
#   Molecule specification                                                 #
#                                                                          #
#   ($Optional)                                                            #
#                                                                          #
#   ($Freeze)                                                              #
#                                                                          #
#   (radii = $Radii)                                                       #
#                                                                          #
############################################################################
        
        def __init__(self,file,Ginput,MolSpec):
                print "\nWriting", file+"_"+Ginput.Append+".mop\n"
                fileout = open(file+"_"+Ginput.Append+".mop", "w")
                fileout.write(""+Ginput.Route+"\n")
                fileout.write(Ginput.Title+"\n\n")
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
  

#Create Orca inp file
class writeOinput:

############################################################################
#                                OrcaG Format                              #
#   %chk = $Link0                                                          #
#   (%mem = $Mem)                                                          #
#   (%nproc = $Nproc)                                                      #
#   Route - desired calculation type                                       #
#                                                                          #
#   $Charge $Mult                                                          #
#   Molecule specification                                                 #
#                                                                          #
############################################################################

        def __init__(self,file,Ginput,MolSpec):
                print "\nWriting", file+"_"+Ginput.Append+".inp\n"
                fileout = open(file+"_"+Ginput.Append+".inp", "w")
                #fileout.write("%chk="+Ginput.Link0+"\n")
                if hasattr(Ginput, "Mem"): fileout.write("%mem="+Ginput.Mem+"\n")
                if hasattr(Ginput, "Nproc"): fileout.write("%nprocshared="+Ginput.Nproc+"\n")
                if hasattr(Ginput, "Linda"): fileout.write("%nproclinda="+Ginput.Linda+"\n")
                fileout.write("# "+Ginput.Title+"\n")
		fileout.write("! "+Ginput.Route+"\n")
                fileout.write("* xyz "+str(MolSpec.CHARGE)+" "+str(MolSpec.MULT)+"\n")
                for i in range(0,MolSpec.NATOMS):
                        fileout.write(MolSpec.ATOMTYPES[i])
                        for j in range(0,3):
                                fileout.write("  "+str(Decimal(str((MolSpec.CARTESIANS[i][j])))))
                        fileout.write("\n")
                fileout.write("*\n")
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


		
