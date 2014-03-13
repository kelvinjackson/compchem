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
#                        weightedNMR.py                       #
#                                                             #
#                      For NICS zz values                     #
#                                                             #
###############################################################

#Python Libraries
import subprocess, sys, os, math


#Chemistry Libaries
from ChemUtils import *
from numpy import *
chemshiftatom=[]
chemshiftatom=raw_input("Chemical shift of which atom(s) number? (Press ENTER for coordinates of ghosts)\n")
if chemshiftatom != "":
        chemshiftatom=float(chemshiftatom)
        chemshiftatom2=raw_input("to...?\n")
else: chemshiftatom2 =""
if chemshiftatom2 != "":
    chemshiftatom2=float(chemshiftatom2)
#nicsplus=raw_input("Plot NICS perpendicular to plane?")
#ringsize=int(float(ringsize))
#ringsize = ringsize-1
#print ringsize
#Read molecule data from an output file


class getoutData:

        def __init__(self, file):

                if not os.path.exists(file+".out"):
                        if not os.path.exists(file+".log"):
                                print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)


                def getFORMAT(self, outlines):
                        for i in range(0,len(outlines)):
                                if outlines[i].find("Gaussian") > -1: self.FORMAT = "Gaussian";
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
                        anharmonic_geom=0
                        outtest =0
                        outtest2=0
                        for i in range(0,len(outlines)):
                            if outlines[i].find("Standard orientation") > -1:
                                standor = i
                                outtest2 = 1
                            #if outtest2 == 1:
                            if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1 and outtest2 == 1:
                                self.NATOMS = i-standor-6
                                outtest = 1
                            #if outtest == 0:
                            if outlines[i].find("Input orientation") > -1 and outtest == 0:
                                standor2 = i
                            if outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1 and outtest == 0:
                                self.NATOMS = i-standor2-6
                            #if outtest==0:
                            #if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1:
                                #self.NATOMS = i-standor2-6
                        try: standor
                        except NameError: pass
                        else:
                            for i in range (standor+5,standor+5+self.NATOMS):
                                self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
                                self.CARTESIANS.append([float(outlines[i].split()[3]),float(outlines[i].split()[4]),float(outlines[i].split()[5])])
                        try: standor2
                        except NameError: pass
                        else:
                            for i in range (standor2+5,standor2+5+self.NATOMS):
                                self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
                                self.CARTESIANS.append([float(outlines[i].split()[3]),float(outlines[i].split()[4]),float(outlines[i].split()[5])])
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

def find_centroid(ringatoms,fileData):
    xtot = 0
    xvals=[]
    yvals=[]
    zvals=[]
    for x in ringatoms:
        #print fileData.CARTESIANS[x][0]
        xtot = xtot + fileData.CARTESIANS[x][0]
        xvals.append(fileData.CARTESIANS[x][0])
        #print xvals
    xav = xtot/ringsize
    ytot = 0
    for x in ringatoms:
        #print fileData.CARTESIANS[x][0]
        ytot = ytot + fileData.CARTESIANS[x][1]
        yvals.append(fileData.CARTESIANS[x][1])
    yav = ytot/ringsize
    ztot = 0
    for x in ringatoms:
        #print fileData.CARTESIANS[x][0]
        ztot = ztot + fileData.CARTESIANS[x][2]
        zvals.append(fileData.CARTESIANS[x][2])
    zav = ztot/ringsize
        #x_av = fileData.CARTESIANS[atomA][0]
    print "Centroid at:", xav, yav, zav  #gives position of centroid
    #print xvals, yvals, zvals
    return xvals, yvals, zvals, xav, yav, zav
    
def find_coeffplane(ringatoms, fileData):
    rotated = 0
    xvals, yvals, zvals, xav, yav, zav = find_centroid(ringatoms, fileData)
    print xvals
    print yvals
    print zvals
    xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum = get_squares_list(ringatoms, xvals, yvals, zvals)
    if xsum == 0.0 and ysum == 0.0:
        rotated = 3
        print "Can't define a ring by points in a line"
        print "This is going to go horribly wrong"
    if xsum == 0.0:
        new_xvals = yvals
        new_yvals = zvals
        new_zvals = xvals
        xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum = get_squares_list(ringatoms, xvals, yvals, zvals)
        rotated = 1
    if ysum == 0.0:
        new_xvals = zvals
        new_yvals = xvals
        new_zvals = yvals
        xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum = get_squares_list(ringatoms, xvals, yvals, zvals)
        rotated = 2
        
    coeffplane = do_matrix_stuff(xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum, ringatoms)
    
    return coeffplane, xav, yav, zav, rotated
        
        
        #print "\nSuccessfully read geometry output", file
                
def get_squares_list(ringatoms, xvals, yvals, zvals):
####################Necessary summations
    xysum = 0
    y2sum = 0
    x2sum = 0
    zsum = 0
    ysum = 0
    xsum = 0
    xzsum = 0
    yzsum = 0
    for n in range(len(ringatoms)):
        xy = xvals[n]*yvals[n]
        #print xy
        xysum = xy+xysum
        #print xysum
        xz = xvals[n]*zvals[n]
        #print xy
        xzsum = xz+xzsum
        yz = yvals[n]*zvals[n]
            #print xy
        yzsum = yz+yzsum
        #print xzsum
        x = xvals[n]
        #print xy
        xsum = x+xsum
        #print xsum
        y = yvals[n]
        #print xy
        ysum = y+ysum
        #print ysum
        z = zvals[n]
        #print xy
        zsum = z+zsum
        #print zsum
        x2 = xvals[n]*xvals[n]
        #print xy
        x2sum = x2+x2sum
        #print x2sum
        y2 = yvals[n]*yvals[n]
        #print xy
        y2sum = y2+y2sum
    return xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum

def do_matrix_stuff(xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum, ringatoms):
###################Matrix and vector used for least squares best fit plane
    a=matrix([[x2sum, xysum, xsum],[xysum, y2sum, ysum],[xsum, ysum, len(ringatoms)]]) #3x3 matrix
    b=matrix([[xzsum],[yzsum],[zsum]]) #3x1 matrix
    #print a
    #print b
    try: coeffplane=a.I*b
    except linalg.linalg.LinAlgError: coeffplane = matrix([[0.0],[0.0],[0.0]]) 
    
    return coeffplane
                
                
if __name__ == "__main__":

        #0 file(s)
        files = []

        # Takes arguments: (1) input file(s)
        if len(sys.argv) > 1:
                for i in range(1,2):
                        files.append(sys.argv[i].split(".")[0])
        else:
                print "\nWrong number of arguments used. Correct format: weightedNMR file(s) ringatom numbers\n"
                sys.exit()
        for file in files:
                ringatoms = []
                ringsize=len(sys.argv)-2
                print "Ring Size =",ringsize
                fileData = getoutData(file)
                #print fileData, file, files
                atomA=int(sys.argv[2])-1
                atomB=int(sys.argv[3])-1
                if ringsize>=3:atomC=int(sys.argv[4])-1	#minimum size for a ring
                if ringsize>=4:atomD=int(sys.argv[5])-1
                if ringsize>=5:atomE=int(sys.argv[6])-1
                if ringsize>=6:atomF=int(sys.argv[7])-1
                if ringsize>=7:atomG=int(sys.argv[8])-1
                if ringsize>=8:atomH=int(sys.argv[9])-1
                if ringsize>=9:atomI=int(sys.argv[10])-1
                if ringsize>=10:atomJ=int(sys.argv[11])-1
                ringatoms.append(atomA)#;print "atom A", fileData.CARTESIANS[atomA]
                ringatoms.append(atomB)#;print "atom B", fileData.CARTESIANS[atomB]
                if ringsize>=3:ringatoms.append(atomC)#;print "atom C", fileData.CARTESIANS[atomC]
                if ringsize>=4:ringatoms.append(atomD)#;print "atom D", fileData.CARTESIANS[atomD]
                if ringsize>=5:ringatoms.append(atomE)#;print "atom E", fileData.CARTESIANS[atomE]
                if ringsize>=6:ringatoms.append(atomF)#;print "atom F", fileData.CARTESIANS[atomF]
                if ringsize>=7:ringatoms.append(atomG)#;print "atom G", fileData.CARTESIANS[atomG]
                if ringsize>=8:ringatoms.append(atomH)#;print "atom H", fileData.CARTESIANS[atomH]
                if ringsize>=9:ringatoms.append(atomI)#;print "atom I", fileData.CARTESIANS[atomG]
                if ringsize>=10:ringatoms.append(atomJ)#;print "atom J", fileData.CARTESIANS[atomH]
                #ringatoms.append(atomA)
                #print ringatoms
                coeffplane, xav, yav, zav, rotated = find_coeffplane(ringatoms,fileData)
				#horrible guess
                
                #coeffplane=a.I*b #Since ax=b, and x conatains coefficients for best fit plane in for z=ax+by+
                #print '\n'.join(''.join(str(cell) for cell in row) for row in a)
                #print coeffplane
                xcoeff= coeffplane.tolist()[0][0]
                ycoeff= coeffplane.tolist()[1][0]
                cval= coeffplane.tolist()[2][0]

                print "Equation of best-fit plane:","z="+str(xcoeff)+"x+"+str(ycoeff)+"y+"+str(cval)	#This gives the equation for the plane of best-fit
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
                print "Unit vector:", x, y, z #The length of this vector is 1
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
                    
                if not os.path.exists(file+".out"):
                    if not os.path.exists(file+".log"):
                        print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)
                if os.path.exists(file+".out"):
                    outfile = open(file+".out","r") 
                else:
                    outfile = open(file+".log","r")
                outlines = outfile.readlines()
                NMR = []
                NATOMS = len(fileData.CARTESIANS)
                #print chemshiftatom
                if chemshiftatom != "":
                    for i in range(0,len(outlines)):
                        if outlines[i].find(" SCF GIAO Magnetic shielding tensor (ppm):") > -1:
                            j = 0
                            dist=-4 #######Need to change manually - will fix at some point...
                            while j < NATOMS*5:
                                item = {}
                                item['atom_index'] = int(outlines[i+1+j].split()[0])
                                item['elementID'] = outlines[i+1+j].split()[1]
                                item['isotropic'] = float(outlines[i+1+j].split()[4])
                                item['anisotropy'] = float(outlines[i+1+j].split()[7])
                                item['xx'] = float(outlines[i+2+j].split()[1])
                                item['yx'] = float(outlines[i+2+j].split()[3])
                                item['zx'] = float(outlines[i+2+j].split()[5])
                                item['xy'] = float(outlines[i+3+j].split()[1])
                                item['yy'] = float(outlines[i+3+j].split()[3])
                                item['zy'] = float(outlines[i+3+j].split()[5])
                                item['xz'] = float(outlines[i+4+j].split()[1])
                                item['yz'] = float(outlines[i+4+j].split()[3])
                                item['zz'] = float(outlines[i+4+j].split()[5])
                                item['eigenvalues'] = [float(outlines[i+5+j].split()[1]), float(outlines[i+5+j].split()[2]), float(outlines[i+5+j].split()[3])]
                                NMR.append(item)
                                j += 5

                                if chemshiftatom2 == "":
                                    if j/5 == chemshiftatom:
                                        NICS = item['xx']*x+item['yy']*y+item['zz']*z
                                        print "NICS 0 =",NICS
                                else:
                                    for q in range(int(chemshiftatom), int(chemshiftatom2)+1):
                                        if j/5==q:
                                            NICS = item['xx']*x+item['yy']*y+item['zz']*z
                                            print "NICS atom"+str(q),dist,NICS
                                            dist=dist+0.25

                     #                       if j/5==chemshiftatom:


##################NICS trajectory above and below plane
                #if nicsplus == "y" or "Y" or "yes":
                if chemshiftatom=="":
                    w=0
                    range_angstrom=4 ###########################Change this value to adjust range of ghost atom points
                    range_number=range_angstrom*8
                    for w in range(0,range_number+1):
                        v=w-(range_number/2)
                        scalefactor=float(v)/4;# print scalefactor
                        xscale=x*scalefactor
                        yscale=y*scalefactor
                        zscale=z*scalefactor
                        print "Bq", xav+xscale, yav+yscale, zav+zscale






















