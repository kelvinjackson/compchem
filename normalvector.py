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
###############################################################

#Python Libraries
import subprocess, sys, os, math

#Chemistry Libaries
from ChemUtils import *
from numpy import *

## Gets data from Gaussian formatted output file
class getoutData:
	def __init__(self, file):
		
		if not os.path.exists(file+".out"):
			if not os.path.exists(file+".log"):
				print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)
		
		
		def getFORMAT(self, outlines):
			for i in range(0,len(outlines)):
				if outlines[i].find("Gaussian") > -1: self.FORMAT = "Gaussian"; break
				if outlines[i].find("MOPAC") > -1: self.FORMAT = "Mopac"; break
				
		def getATOMTYPES(self, outlines, format):
			self.ATOMTYPES = []
			self.CARTESIANS = []
			anharmonic_geom=0
			outtest =0; outtest2=0
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
					
		if os.path.exists(file+".out"): outfile = open(file+".out","r")
		else: outfile = open(file+".log","r")
		outlines = outfile.readlines()
		getFORMAT(self, outlines)
		getATOMTYPES(self, outlines, self.FORMAT)


def find_centroid(ringatoms,fileData):
	xtot = 0; xvals=[]; yvals=[]; zvals=[]
	for x in ringatoms:
		print fileData.CARTESIANS[x]
		xtot = xtot + fileData.CARTESIANS[x][0]
		xvals.append(fileData.CARTESIANS[x][0])
	xav = xtot/ringsize
	ytot = 0
	for x in ringatoms:
		ytot = ytot + fileData.CARTESIANS[x][1]
		yvals.append(fileData.CARTESIANS[x][1])
	yav = ytot/ringsize
	ztot = 0
	for x in ringatoms:
		ztot = ztot + fileData.CARTESIANS[x][2]
		zvals.append(fileData.CARTESIANS[x][2])
	zav = ztot/ringsize

	print "Centroid at:", xav, yav, zav  #gives position of centroid
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
	# Takes arguments: (1) Gaussian output file(s)
	files = []
	if len(sys.argv) > 1:
		for i in range(1,2): files.append(sys.argv[i].split(".")[0])
	else:
		print "\nWrong number of arguments used. Correct format: weightedNMR file(s) ringatom numbers\n"; sys.exit()

	for file in files:
		ringatoms = []; ringsize=len(sys.argv)-2
		print "Ring Size =",ringsize

		## Get coordinates
		fileData = getoutData(file)

		## Establish the array of atoms in the ring
		for atom in sys.argv[2:]: ringatoms.append(int(atom)-1)
		coeffplane, xav, yav, zav, rotated = find_coeffplane(ringatoms,fileData)
		xcoeff= coeffplane.tolist()[0][0]; ycoeff= coeffplane.tolist()[1][0]; cval= coeffplane.tolist()[2][0]

		print "Equation of best-fit plane:","z="+str(xcoeff)+"x+"+str(ycoeff)+"y+"+str(cval)	#This gives the equation for the plane of best-fit
		####################Make unit vector
		rawvector=array([xcoeff,ycoeff,-1]) #Need to make into unit vector
		x=float(rawvector[0]); y=float(rawvector[1]); z=float(rawvector[2])
		normfactor=1/(x**2+y**2+z**2)**0.5
		x=x*normfactor; y=y*normfactor; z=z*normfactor
		if z<0: z=-z;y=-y;x=-x #Sign flip if z is negative
		print "Unit vector:", x, y, z #The length of this vector is 1

		if rotated == 1:
			print "************ coordinated system was rotated! ***********"
			old_x = z; old_y = x; old_z = y
			if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
			print "Unit vector:", old_x, old_y, old_z
			x = old_x; y = old_y; z = old_z
		if rotated == 2:
			print "************ coordinated system was rotated! ***********"
			old_x = y; old_y = z; old_z = x
			if old_z<0: old_z=-old_z;old_y=-old_y;old_x=-old_x
			print "Unit vector:", old_x, old_y, old_z
			x = old_x; y = old_y; z = old_z
		if rotated == 3:
			print "didn't I tell you this was a bad idea?"
		
