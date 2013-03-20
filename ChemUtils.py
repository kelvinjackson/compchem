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
#                       ChemUtils.py                          #
#   elementID - returns element symbol from atomic number     #
#   atomicNum - returns atomic number from elemental symbol   #
#   bondiRadius - return Bondi radius from atomic number      #
#   digitalMonth - returns month as number 1-12 from name     #
#   reName - replaces occurences of oldname to newname for    #
#            all files in a directory                         #
###############################################################


#Python libararies
import sys, os, commands, string, math


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
		
