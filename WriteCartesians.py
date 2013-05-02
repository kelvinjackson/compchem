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
#                                                             #
#                         ccParse.py                          #
#                                                             #
#                Reads compchem job file(s)                   #
#                                                             #
###############################################################

#Python Libraries 
import subprocess, sys, os

if __name__ == "__main__":
        
        #0 file(s)
        files = []

        # Takes arguments: (1) input file(s)
        if len(sys.argv) > 1: 
                for i in range(1,len(sys.argv)):
                        files.append(sys.argv[i])
        else:
                print "\nWrong number of arguments used. Correct format: ccParse file(s)\n"
                sys.exit()

	for file in files:
		#print "WriteCartesians", file
		print "ccParse.py", file, ">", str(file+"_energy1")
		print "GaussianPrep.py", file, "-route Cartesians -append cart > /dev/null"
		print "cat", file+"_energy1"
		print "rm", file+"_energy1"
		#print "cat", file+"_energy2"
		#print "cat", file+"_energy3"
		print "tail -n +7", file.split(".")[0]+"_cart.com"
		print "rm", file.split(".")[0]+"_cart.com"

