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
#                        PICKNICS                             #
#                                                             #
#    Resolves NICS values from g09 calculations into their    #
#    in-plane and out-of-plane components. For a (specified)  #
#    ring of atoms the least-squares plane and normal vector  #
#    are calculated. If NCS values have been calculated this  #
#    is also performed for selected NBOs.                     #
#                                                             #
#    Alex D'Anthony, Kelvin Jackson and RS Paton, April 2014  #
###############################################################

#Python Libraries
import subprocess, sys, os, math
from numpy import *

periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
        "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

def elementID(massno):
                if massno < len(periodictable): return periodictable[massno]    
                else: return "XX"

class getoutData:
        def __init__(self, file):

		if not os.path.exists(file+".out"):
                        if not os.path.exists(file+".log"):
                                print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)

                def getFORMAT(self, outlines):
                        for i in range(0,len(outlines)):
                                if outlines[i].find("Gaussian") > -1: self.FORMAT = "Gaussian";

                def getTERMINATION(self, outlines,format):
                        if format == "Gaussian":
                                for i in range(0,len(outlines)):
                                        if outlines[i].find("Normal termination") > -1:
                                                self.TERMINATION = "normal"

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
                            if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1 and outtest2 == 1:
                                self.NATOMS = i-standor-6
                                outtest = 1
                            if outlines[i].find("Input orientation") > -1 and outtest == 0:
                                standor2 = i
                            if outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1 and outtest == 0:
                                self.NATOMS = i-standor2-6
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
                
		def getSHIELDING(self,outlines,format):
			self.NMR = []
                    	for i in range(0,len(outlines)):
                        	if outlines[i].find(" SCF GIAO Magnetic shielding tensor (ppm):") > -1:
                            		j = 0
                            		dist=-4 #######Need to change manually - will fix at some point...
                            		while j < self.NATOMS*5:
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
                                		self.NMR.append(item)
                                		j += 5

		def checkNBOS(self,outlines,format):
			nbos = "off"
			for i in range(0,len(outlines)):
				if outlines[i].find("NATURAL BOND ORBITAL ANALYSIS:") > -1: 
					nbos = "on"
			return nbos

                def getNBOS(self,outlines,format):
			self.NBOS = []
			self.NBOtype = []
                        for i in range(0,len(outlines)):
                                if outlines[i].find("------------------ Lewis ----------------------------") > -1: startnbos = i
				if outlines[i].find("---------------- non-Lewis ----------------------------------------------------") > -1: endnbos = i
			for i in range(startnbos, endnbos):
				occupiednbos = ["CR", "LP", "BD"]
				for occ in occupiednbos:
					if outlines[i].find(occ) > -1:
						self.NBOS.append(outlines[i][23:36].lstrip().rstrip()) 
						self.NBOtype.append(outlines[i][16:23].lstrip().rstrip())	

		def getNCS(self,outlines,format):
                        self.NCS = []
			n = -1
                        for i in range(0,len(outlines)):
                                if outlines[i].find("Full Cartesian NMR shielding tensor (ppm) for atom") > -1:
					atomtype = outlines[i].split()[8].split("(")[0]
				if outlines[i].find("L  NL    XX      XY      XZ      YX      YY      YZ      ZX      ZY      ZZ") > -1:
					for j in range(0,200):
						if outlines[i+2+j].find("Total") > -1:
							endline = j-4; break
				
					self.NCS.append([])
					n = n + 1
					for j in range(0,endline):
						item = {}
						item['elementID'] =  atomtype
						#print outlines[i+2+j][:5]
						if len(outlines[i+2+j][:5].split()) > 0:
                                                	item['nbo'] = int(float(outlines[i+2+j].split()[0]))
						else: item['nbo'] = 0
						item['xx'] = float(outlines[i+2+j].split()[1])
                                                item['xy'] = float(outlines[i+2+j].split()[2])
                                                item['xz'] = float(outlines[i+2+j].split()[3])
                                                item['yx'] = float(outlines[i+2+j].split()[4])
						item['yy'] = float(outlines[i+2+j].split()[5])
                                                item['yz'] = float(outlines[i+2+j].split()[6])
                                                item['zx'] = float(outlines[i+2+j].split()[7])
                                                item['zy'] = float(outlines[i+2+j].split()[8])
                                                item['zz'] = float(outlines[i+2+j].split()[9])
						
						self.NCS[n].append(item)

		if os.path.exists(file+".out"):outfile = open(file+".out","r")
                else: outfile = open(file+".log","r")
                outlines = outfile.readlines()

                getFORMAT(self, outlines)
                getTERMINATION(self, outlines,self.FORMAT)
                if self.TERMINATION == "normal":
                	getATOMTYPES(self, outlines, self.FORMAT)
			getSHIELDING(self,outlines,self.FORMAT)
			
		if checkNBOS(self, outlines, self.FORMAT) == "on": 
			getNBOS(self,outlines,self.FORMAT)
			#print self.NBOS
			#print self.NBOtype
			getNCS(self,outlines,self.FORMAT)
	
def find_centroid(ringatoms,fileData):
        xtot = 0; xvals=[]; yvals=[]; zvals=[]
        for x in ringatoms:
                #print "CARTS", fileData.CARTESIANS[x]
                xtot = xtot + fileData.CARTESIANS[x][0]
                xvals.append(fileData.CARTESIANS[x][0])
        xav = xtot/len(ringatoms)
        ytot = 0
        for x in ringatoms:
                ytot = ytot + fileData.CARTESIANS[x][1]
                yvals.append(fileData.CARTESIANS[x][1])
        yav = ytot/len(ringatoms)
        ztot = 0
        for x in ringatoms:
                ztot = ztot + fileData.CARTESIANS[x][2]
                zvals.append(fileData.CARTESIANS[x][2])
        zav = ztot/len(ringatoms)

        #print "Centroid at:", xav, yav, zav  #gives position of centroid
        return xvals, yvals, zvals, xav, yav, zav


def find_coeffplane(ringatoms, fileData):
        rotated = 0
        xvals, yvals, zvals, xav, yav, zav = find_centroid(ringatoms, fileData)
        #print xvals, yvals, zvals
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


def get_squares_list(ringatoms, xvals, yvals, zvals):
####################Necessary summations
        xysum = 0; y2sum = 0; x2sum = 0; zsum = 0; ysum = 0; xsum = 0; xzsum = 0; yzsum = 0
        for n in range(len(ringatoms)):
                xy = xvals[n]*yvals[n]
                xysum = xy+xysum
                xz = xvals[n]*zvals[n]
                xzsum = xz+xzsum
                yz = yvals[n]*zvals[n]
                yzsum = yz+yzsum
                x = xvals[n]
                xsum = x+xsum
                y = yvals[n]
                ysum = y+ysum
                z = zvals[n]
                zsum = z+zsum
                x2 = xvals[n]*xvals[n]
                x2sum = x2+x2sum
                y2 = yvals[n]*yvals[n]
                y2sum = y2+y2sum
        return xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum

def do_matrix_stuff(xzsum, xysum, xsum, ysum, zsum, x2sum, y2sum, yzsum, ringatoms):
        ###################Matrix and vector used for least squares best fit plane
        a=matrix([[x2sum, xysum, xsum],[xysum, y2sum, ysum],[xsum, ysum, len(ringatoms)]]) #3x3 matrix
        b=matrix([[xzsum],[yzsum],[zsum]]) #3x1 matrix
        try: coeffplane=a.I*b
        except linalg.linalg.LinAlgError: coeffplane = matrix([[0.0],[0.0],[0.0]])
        return coeffplane


if __name__ == "__main__":
	distlist = []
	files = []
        # Takes arguments: (1) input file(s)
        if len(sys.argv) > 1:
                for i in range(1,2):
                        files.append(sys.argv[i].split(".")[0])
        else:
                print "\nWrong number of arguments used. Correct format: weightedNMR file(s) ringatom numbers\n"
                sys.exit()

        for file in files:
                print "######################################################################"
                ringatoms = []; ringsize=len(sys.argv)-2
                #print "Ring Size =",ringsize

                ## Get coordinates
                fileData = getoutData(file)

                ## Establish the array of atoms in the ring
                for atom in sys.argv[2:]: ringatoms.append(int(atom)-1)
                coeffplane, xav, yav, zav, rotated = find_coeffplane(ringatoms,fileData)
                xcoeff= coeffplane.tolist()[0][0]; ycoeff= coeffplane.tolist()[1][0]; cval= coeffplane.tolist()[2][0]

                #print "Equation of best-fit plane:","z="+str(xcoeff)+"x+"+str(ycoeff)+"y+"+str(cval)    #This gives the equation for the plane of best-fit
                ####################Make unit vector
                rawvector=array([xcoeff,ycoeff,-1]) #Need to make into unit vector
                x=float(rawvector[0]); y=float(rawvector[1]); z=float(rawvector[2])
                normfactor=1/(x**2+y**2+z**2)**0.5
                x=x*normfactor; y=y*normfactor; z=z*normfactor
                if z<0: z=-z;y=-y;x=-x #Sign flip if z is negative
                #print "Unit vector:", x, y, z #The length of this vector is 1
        
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
                   
                for tensor in fileData.NMR:
			
			if tensor['elementID'] == "Bq": 
				xdist = fileData.CARTESIANS[(tensor['atom_index']-1)][0] - xav
				ydist = fileData.CARTESIANS[(tensor['atom_index']-1)][1] - yav
				zdist = fileData.CARTESIANS[(tensor['atom_index']-1)][2] - zav
				dist = xdist*x + ydist*y + zdist*z
				distlist.append(dist)
				print "NICS @",
				print ("%.2f" % dist+"   iso:").rjust(13),
				nics_iso = tensor['isotropic'] * -1
				nics_oop = (tensor['xx']*x+tensor['yy']*y+tensor['zz']*z)* -1
				nics_ip = (nics_iso * 3 - nics_oop) / 2
				print str(nics_iso).rjust(7),
				print "  O-O-P (zz):", ("%.3f" % nics_oop).rjust(7) ,"  I-P (zz):", ("%.3f" % nics_ip).rjust(7)
		print "######################################################################"

		if hasattr(fileData, "NBOS"):
			orblist = []
			NBOLIST =raw_input("Natural bonding orbital number(s) of interest...? ")                 
   			for orb in NBOLIST.split(): orblist.append(int(orb))          
			if len(orblist) == 0: orblist = list(xrange(len(fileData.NBOS)))
			print "######################################################################"
			n = 0
			print "NBOs:          ",
			for ORB in orblist: print (fileData.NBOS[ORB-1].replace(' ', '')+":"+fileData.NBOtype[ORB-1].split()[0]).ljust(15),
			print "SUM".ljust(15),
			combined = 0.0
			for naturaltensor in fileData.NCS:
				if len(naturaltensor) > 0:
					if naturaltensor[0]['elementID'] == "gh":
						combined = 0.0
						print "\nNICS @",
						print ("%.2f" % distlist[n]).rjust(6),
						n = n+1
						for ORB in orblist:
							present="absent"
							for i in range(0,len(naturaltensor)):
								if naturaltensor[i]['nbo'] == ORB: 
									present="correct"
									#print naturaltensor[i]['xx']*-0.333+naturaltensor[i]['yy']*-0.333+naturaltensor[i]['zz']*-0.333,
									combined = combined + (naturaltensor[i]['xx']*-x+naturaltensor[i]['yy']*-y+naturaltensor[i]['zz']*-z)
									print "  MOzz:", ("%.3f" % (naturaltensor[i]['xx']*-x+naturaltensor[i]['yy']*-y+naturaltensor[i]['zz']*-z)).rjust(7),
							if present != "correct":print "  MOzz:", ("%.3f" % 0.00).rjust(7),
						print  ("%.3f" % combined).rjust(7),
			print "\n######################################################################"
