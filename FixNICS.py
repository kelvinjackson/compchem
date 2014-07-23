#!/usr/bin/python

from GaussianPrep import *
import sys
from PickNICS import *
import ccParse

if __name__ == "__main__":
        #0 file(s)
        file = sys.argv[1].split(".")[0]
        job =[]
        # Takes arguments: (1) input file(s)
        if len(sys.argv) > 1:
            for i in range(1,len(sys.argv)):
                if sys.argv[i][0:1] == "-" and sys.argv[i][0:3] != "--L":
                    job.append([sys.argv[i],sys.argv[i+1]])
        else:
            print "\nWrong number of arguments used. Correct format: \n"
            sys.exit()

        if os.path.exists(file+".out") or os.path.exists(file+".log"): fileData = ccParse.getoutData(file)
        if os.path.exists(file+".com"): fileData = ccParse.getinData(file)
        
	NICS_RANGE = 4
        for item in job:
            if "NICSRANGE" in item[0].upper():
                NICS_RANGE = int(item[1])
        NICS_NUMBER = 16
        for item in job:
            if "NICSNUMBER" in item[0].upper():
                NICS_NUMBER = int(item[1])
        for item in job:
            if "RING" in item[0].upper():
                ringatoms = []
                for atom in item[1].split():
                    ringatoms.append(int(atom)-1)
                #print ringatoms
                coeffplane, xav, yav, zav, rotated = find_coeffplane(ringatoms,fileData)
                #horrible guess
                
                #coeffplane=a.I*b #Since ax=b, and x conatains coefficients for best fit plane in for z=ax+by+
                #print '\n'.join(''.join(str(cell) for cell in row) for row in a)
                #print coeffplane
                    
                xcoeff= coeffplane.tolist()[0][0]
                ycoeff= coeffplane.tolist()[1][0]
                cval= coeffplane.tolist()[2][0]

#print "Equation of best-fit plane:","z="+str(xcoeff)+"x+"+str(ycoeff)+"y+"+str(cval)	#This gives the equation for the plane of best-fit
####################Make unit vector
                rawvector=array([xcoeff,ycoeff,-1]) #Need to make into unit vector
               # print rawvector
                x=float(rawvector[0])
                y=float(rawvector[1])
                z=float(rawvector[2])
                #print x,y,z
                normfactor=1/(x**2+y**2+z**2)**0.5
                x=x*normfactor; y=y*normfactor; z=z*normfactor
                if z<0: z=-z;y=-y;x=-x #Sign flip if z is negative
                #print "Unit vector:", x, y, z #The length of this vector is 1
                #print "NICS 1 point", x+xav, y+yav, z+zav
                if rotated == 1:
                    print "************ coordinated system was rotated! ***********"
                    old_x = z
                    old_y = x
                    old_z = y
                    if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
                    print "Unit vector:", old_x, old_y, old_z
                    x = old_x
                    y = old_y
                    z = old_z
                if rotated == 2:
                    print "************ coordinated system was rotated! ***********"
                    old_x = y
                    old_y = z
                    old_z = x
                    if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
                    print "Unit vector:", old_x, old_y, old_z
                    x = old_x
                    y = old_y
                    z = old_z
                if rotated == 3:
                    print "didn't I tell you this was a bad idea?"

                spacing = float(NICS_RANGE)/float(NICS_NUMBER)
                for w in range(-NICS_NUMBER,NICS_NUMBER+1):
                    scalefactor = w*spacing
                    xscale=x*scalefactor
                    yscale=y*scalefactor
                    zscale=z*scalefactor
                    fileData.NATOMS += 1
                    fileData.ATOMTYPES.append("Bq")
                    fileData.CARTESIANS.append([xav+xscale, yav+yscale, zav+zscale])
                        #print "Bq", xav+xscale, yav+yscale, zav+zscale
        #print fileData.ATOMTYPES
        Ginput = getGinput(file, job)
        Gwrite = writeGinput(file, Ginput, fileData)
