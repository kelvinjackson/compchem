#!/usr/bin/python

import sys
import string

nfiles =len(sys.argv)-2
filename = sys.argv[1]
oldword = sys.argv[nfiles]
newword = sys.argv[nfiles+1]

for i in range(1,nfiles):
        filename=sys.argv[i]
        if filename.find(oldword)>-1:
		newfilename=string.replace(filename, oldword, newword)
		print "mv ",filename, newfilename
	#tempfile = "temp_"+filename
        #print "sed -e \"s/"+oldword+"/"+newword+"/g\"", filename,">",tempfile
        #print "mv ",tempfile, filename

