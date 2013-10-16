#!/usr/bin/python

import sys
import string
import subprocess

nfiles =len(sys.argv)-2
filename = sys.argv[1]
oldword = sys.argv[nfiles]
newword = sys.argv[nfiles+1]
print nfiles

for i in range(1,nfiles):
        filename=sys.argv[i]
        if filename.find(oldword)>-1:
		newfilename=string.replace(filename, oldword, newword)
		command =  "mv "+filename+"  "+newfilename
		try:
                	retcode = subprocess.call(command, shell=True)
                	#if retcode != 0: return -1
                	#else: return 1
        	except OSError, e:
                	print >>sys.stderr, log.Write("\nERROR")
                	#return -1
######
	#tempfile = "temp_"+filename
        #print "sed -e \"s/"+oldword+"/"+newword+"/g\"", filename,">",tempfile
        #print "mv ",tempfile, filename

