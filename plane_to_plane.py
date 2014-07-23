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
from numpy import *
#Some useful arrays
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
	"Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

calendar=["","jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"]

atomicmass = [1.008, 4.003, 6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07, 35.45, 39.95, 39.10, 40.08, 44.96, 47.87, 50.94, 52.00, 54.94, 55.84, 58.93, 58.69, 
	63.55, 65.39, 69.72, 72.61, 74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91, 91.22, 92.91, 95.94, 99.0, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29]


def elementID(massno):
		if massno < len(periodictable): return periodictable[massno]	
		else: return "XX"


def atomicnumber(element):
	atomicno = 0
	for i in range(0,len(periodictable)):
		if element == periodictable[i]: atomicno = i
	return atomicno		


def bondiRadius(massno):
	#Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452, except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391
	#Radii that are not available in either of these publications have RvdW = 2.00 Angstrom
	
	bondi = [0.0,1.09, 1.40, 1.82,2.00,2.00,1.70,1.55,1.52,1.47,1.54,2.27,1.73,2.00,2.10,1.80,1.80,1.75,1.88,2.75,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,1.87,2.00,1.85,1.90,
	1.85,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,2.00,2.06,1.98,2.16,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.72,1.66,1.55,1.96,2.02,2.00,2.00,2.00,
	2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.86]
	if massno<len(bondi): radius = bondi[massno]
	else: radius = 2.0 
	return radius
	

def digitalMonth(month):
	digital = 0
	for i in range(0,len(calendar)):
		if calendar[i] in month.lower(): digital = i
	return digital


def calcdist(atoma,atomb,coords):
	x1=coords[atoma][0]
	y1=coords[atoma][1]
	z1=coords[atoma][2]
	x2=coords[atomb][0]
	y2=coords[atomb][1]
	z2=coords[atomb][2]
	ba = [x1-x2, y1-y2, z1-z2]
	dist = math.sqrt(ba[0]*ba[0]+ba[1]*ba[1]+ba[2]*ba[2])
	return dist
	
	
def calcangle(atoma,atomb,atomc,coords):
	x1=coords[atoma][0]
	y1=coords[atoma][1]
	z1=coords[atoma][2]
	x2=coords[atomb][0]
	y2=coords[atomb][1]
	z2=coords[atomb][2]
	x3=coords[atomc][0]
	y3=coords[atomc][1]
	z3=coords[atomc][2]
	ba = [x1-x2, y1-y2, z1-z2]
	bc = [x3-x2, y3-y2, z3-z2]
	angle = 180.0/math.pi*math.acos((ba[0]*bc[0]+ba[1]*bc[1]+ba[2]*bc[2])/(math.sqrt(ba[0]*ba[0]+ba[1]*ba[1]+ba[2]*ba[2])*math.sqrt(bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]))) 
	return angle
	

def calcdihedral(atoma,atomb,atomc,atomd,coords):
	x1=coords[atoma][0]
	y1=coords[atoma][1]
	z1=coords[atoma][2]
	x2=coords[atomb][0]
	y2=coords[atomb][1]
	z2=coords[atomb][2]
	x3=coords[atomc][0]
	y3=coords[atomc][1]
	z3=coords[atomc][2]
	x4=coords[atomd][0]
	y4=coords[atomd][1]
	z4=coords[atomd][2]
	ax= (y2-y1)*(z2-z3)-(z2-z1)*(y2-y3)
	ay= (z2-z1)*(x2-x3)-(x2-x1)*(z2-z3)
	az= (x2-x1)*(y2-y3)-(y2-y1)*(x2-x3)
	bx= (y3-y2)*(z3-z4)-(z3-z2)*(y3-y4)
	by= (z3-z2)*(x3-x4)-(x3-x2)*(z3-z4)
	bz= (x3-x2)*(y3-y4)-(y3-y2)*(x3-x4)
	nbx= (y2-y3)*(z4-z3)-(z2-z3)*(y4-y3)
	nby= (z2-z3)*(x4-x3)-(x2-x3)*(z4-z3)
	nbz= (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)
	torsion=180.0/math.pi*math.acos((ax*bx+ay*by+az*bz)/(math.sqrt(ax*ax+ay*ay+az*az)*math.sqrt(bx*bx+by*by+bz*bz)))
	sign=180.0/math.pi*math.acos((nbx*(x2-x1)+nby*(y2-y1)+nbz*(z2-z1))/(math.sqrt(nbx*nbx+nby*nby+nbz*nbz)*math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))))
	if sign<90.0:
		torsion=torsion*-1.0
	return torsion
	
	
def reName(dir, oldname, newname):
	files=commands.getstatusoutput("ls "+dir+"/*"+oldname+"*")
	if files[0] == 0:
		for file in string.split(files[1],'\n'):
			newfile = file.replace(oldname, newname)
			print "Renaming",file,"to",newfile
			commands.getoutput("mv "+file+" "+newfile) 
	else:
		print "No files Found"
	
## Gets data from Gaussian formatted output file
class getoutData:
	def __init__(self, file):
		
		if not os.path.exists(file):
			print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)
		
		def getATOMTYPES(self, outlines, format):
			self.ATOMTYPES = []
			self.CARTESIANS = []
			anharmonic_geom=0
			outtest =0; outtest2=0
			self.CARTESIANS.append([0.0,0.0,0.0])
			self.ATOMTYPES.append("null")
			for i in range(0,len(outlines)):
				if outlines[i].find("ATOM") > -1:
					self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
					self.CARTESIANS.append([float(outlines[i][30:38]),float(outlines[i][38:46]),float(outlines[i][46:54])])
					#self.CARTESIANS.append([float(outlines[i].split()[5]),float(outlines[i].split()[6]),float(outlines[i].split()[7])])
					
		if os.path.exists(file):
			outfile = open(file,"r") 	
			outlines = outfile.readlines()
		self.FORMAT = "PDB"
		getATOMTYPES(self, outlines, self.FORMAT)

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
	# Takes arguments: (1) Gaussian output file(s)
	files = []
	if len(sys.argv) > 1:
		for i in range(1,len(sys.argv)): files.append(sys.argv[i])
	else:
		print "\nWrong number of arguments used. Correct format: weightedNMR file(s) ringatom numbers\n"; sys.exit()

	for file in files:
		ringatoms = [1975, 1976, 1977, 1979, 1985, 1987]
		ringsize=len(ringatoms)
		#print "Ring Size =",ringsize
		argatoms = [1552,1554,1555,1558]
			
		## Get coordinates
		fileData = getoutData(file)

		## Establish the array of atoms in the ring
		print "AROMATIC RING", ringatoms
		print "ARG ATOMS", argatoms
		coeffplane, xav, yav, zav, rotated = find_coeffplane(ringatoms,fileData)
		xcoeff= coeffplane.tolist()[0][0]; ycoeff= coeffplane.tolist()[1][0]; cval= coeffplane.tolist()[2][0]
		argplane, argxav, argyav,argzav, argrotated = find_coeffplane(argatoms,fileData)
		argxcoeff= argplane.tolist()[0][0]; argycoeff= argplane.tolist()[1][0]; argcval= argplane.tolist()[2][0]
		d_vec = [xav-argxav,yav-argyav,zav-argzav]
		dist = (d_vec[0] ** 2) + (d_vec[1] ** 2) + (d_vec[2] ** 2)
		dist = dist ** 0.5		

		#print "Equation of best-fit plane specfied:","z="+str(xcoeff)+"x+"+str(ycoeff)+"y+"+str(cval)	#This gives the equation for the plane of best-fit
		#print "Equation of arginine best-fit plane:","z="+str(argxcoeff)+"x+"+str(argycoeff)+"y+"+str(argcval)   #This gives the equation for the plane of best-fit
	
		####################Make unit vector
		rawvector=array([xcoeff,ycoeff,-1]) #Need to make into unit vector
		x=float(rawvector[0]); y=float(rawvector[1]); z=float(rawvector[2])
		normfactor=1/(x**2+y**2+z**2)**0.5
		x=x*normfactor; y=y*normfactor; z=z*normfactor
		if z<0: z=-z;y=-y;x=-x #Sign flip if z is negative
		#print "Aromatic Unit vector:", x, y, z #The length of this vector is 1

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
		
		rawargvector=array([argxcoeff,argycoeff,-1]) #Need to make into unit vector
                argx=float(rawargvector[0]); argy=float(rawargvector[1]); argz=float(rawargvector[2])
                argnormfactor=1/(argx**2+argy**2+argz**2)**0.5
                argx=argx*argnormfactor; argy=argy*argnormfactor; argz=argz*argnormfactor
                if argz<0: argz=-argz;argy=-argy;argx=-argx #Sign flip if z is negative
                #print "Arg Unit vector:", argx, argy, argz #The length of this vector is 1
		
		dotprod = x*argx + y*argy + z*argz
		angle_between_planes = 180.0 * math.acos(dotprod) / math.pi
		print file, "ANGLE BETWEEN PLANES (deg)", angle_between_planes
		print file, "DISTANCE BETWEEN COM (Ang)", dist
		print ""
