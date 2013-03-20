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
#                        eval_D3.py                           #
#     Kelvin Jackson & Robert Paton, University of Oxford     #
#                          2013                               #
#         evalute D3 steric repulsion and dispersion          #
###############################################################

# Dependent on parameter file
from pars import *

#Python libararies
import random, sys, os, commands, string, math

## Check for integer when parsing ##
def is_number(s):
    try: int(s); return True
    except ValueError: return False

## Arrays for attractive and repulsive interactions ##
attractive_vdw=[0]
repulsive_vdw=[0]
total_vdw=[0]

## Functional Specific D3 parameters
rs8 = 1.0
repfac= 1.0

## Conversion factors ##
autoang = 0.52917726
autokcal = 627.509541
c6conv=(0.001/2625.4999)/(0.052917726**6)

## Global D3 parameters ##
## Exponents used in distance dependent damping factors for R6, R8 and R10 terms
alpha6 = 14
alpha8 = alpha6 + 2
alpha10 = alpha8 + 2
## Constants used to determine the fractional connectivity between 2 atoms:
## k1 is the exponent used in summation, k2 is used a fraction of the summed single-bond radii
k1 = 16.0
k2 = 4.0/3.0
k3 = -4.0
## D3 is parameterized up to element 94
max_elem = 94
## maximum connectivity
maxc = 5


## Work out the atoms in the same molecules
def getMollist(bondmatrix,startatom):
	
	# The list of atoms in a molecule
	atomlist=[]
	atomlist.append(startatom)

	molecule1=[]
	nextlot=[]
	count = 0

	while count<100:
		nextlot=[]
		for atom in atomlist:
			#print atom, bondmatrix[atom]
			
			for i in range(0,len(bondmatrix[atom])):
				if bondmatrix[atom][i] == 1:
					alreadyfound = 0
					for at in atomlist:
						if i == at: alreadyfound = 1
					if alreadyfound == 0: atomlist.append(i)

		count=count+1	

	return atomlist
	

## DFT derived values for diatomic cutoff radii from Grimme ##
## These are read from pars.py and converted from atomic units into Angstrom
r = [[0]*max_elem for x in xrange(max_elem)]
k=0
for i in range(0,max_elem):
	for j in range(0,i+1):
		r[i][j]=r0ab[k]/autoang
		r[j][i]=r0ab[k]/autoang
		k=k+1

## PBE0/def2-QZVP atomic values for multipole coefficients read from pars.py ##
for i in range(0,max_elem):
	dum=0.5*r2r4[i]*float(i+1)**0.5
	r2r4[i]=math.pow(dum,0.5)


## Reference systems are read in to compute coordination number dependent dispersion coefficients
def copyc6(max_elem, maxc):
	c6ab = [[0]*max_elem for x in xrange(max_elem)]
	nlines = 32385
	
	for iat in range(0,max_elem):
		for jat in range(0,max_elem):
			c6ab[iat][jat]=[[0]*maxc for x in xrange(maxc)]
	
	kk=0
	for nn in range(0,nlines):
		kk=(nn*5)
		iadr=0
		jadr=0
		#print pars[kk], pars[kk+1], pars[kk+2]
		iat=int(pars[kk+1])-1
		jat=int(pars[kk+2])-1
		
		while iat > 99:
			iadr=iadr+1
			iat=iat-100
		while jat > 99:
			jadr=jadr+1
			jat=jat-100
		
		c6ab[iat][jat][iadr][jadr]=[]
		c6ab[iat][jat][iadr][jadr].append(pars[kk])
		c6ab[iat][jat][iadr][jadr].append(pars[kk+3])
		c6ab[iat][jat][iadr][jadr].append(pars[kk+4])
		
		c6ab[jat][iat][jadr][iadr]=[]
		c6ab[jat][iat][jadr][iadr].append(pars[kk])
		c6ab[jat][iat][jadr][iadr].append(pars[kk+4])
		c6ab[jat][iat][jadr][iadr].append(pars[kk+3])
	return c6ab


def getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,a,b):
	for i in range(0,max_elem):
		if atomtype[a].find(elements[i])>-1:iat=i
		if atomtype[b].find(elements[i])>-1:jat=i
	
	c6mem = -1.0E99
	rsum = 0.0
	csum = 0.0
	c6  = 0.0
	for i in range(0,mxc[iat]):
		for j in range(0,mxc[jat]):
			if isinstance(c6ab[iat][jat][i][j], (list, tuple)):
				c6=c6ab[iat][jat][i][j][0]
				if c6>0:
					c6mem=c6
					cn1=c6ab[iat][jat][i][j][1]
					cn2=c6ab[iat][jat][i][j][2]
					
					r=(cn1-cn[a])**2+(cn2-cn[b])**2
					tmp1=math.exp(k3*r)
					rsum=rsum+tmp1
					csum=csum+tmp1*c6
		
	if(rsum>0):
		c6=csum/rsum
	else:
		c6=c6mem	
	return c6

def ncoord(natom, rcov, atomtype, xco, yco, zco, max_elem, autoang, k1, k2):
	cn =[]
	for i in range(0,natom):
		
		xn = 0.0
		for iat in range(0,natom):
			if iat != i:
				dx = xco[iat] - xco[i]
				dy = yco[iat] - yco[i]
				dz = zco[iat] - zco[i]
				r2 = dx*dx+dy*dy+dz*dz
				r = math.pow(r2,0.5)
				
				for k in range(0,max_elem):
					if atomtype[i].find(elements[k])>-1:Zi=k
					if atomtype[iat].find(elements[k])>-1:Ziat=k
				
				rco = rcov[Zi]+rcov[Ziat]
				rco = rco*k2
				rr=rco/r
				damp=1.0/(1.0+math.exp(-k1*(rr-1.0)))
				#print Zi, Ziat,r, rcov[Zi], rcov[Ziat], rco,rr, damp
				xn=xn+damp
		
		cn.append(xn)
	return cn


## Get from pars.py
c6ab = copyc6(max_elem, maxc)


## The periodic table...
elm=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']


## The computation of the classical non-bonding interaction
class calcD3:
	
	def __init__(self, file, s6, rs6, s8):
		
		file = file.split(".com")[0]
		infile = open(file+".com","r")
		inlines = infile.readlines()

		## Get past the standard Gaussian Input lines at the top of the file
		## This may need making more general 
		line0=inlines[0]; line1=inlines[1]; line2=inlines[2]; line3=inlines[3]; line4=inlines[4]; line5=inlines[5]
		for j in range(0,8):
			if len(inlines[j])>0 and len(inlines[j])<130:
				if len(inlines[j].split())==2 or len(inlines[j].split())==6: titleline=j; break

		## Arrays for atoms and Cartesian coordinates ##
		atom=[]
		atomtype=[]
		xco=[]
		yco=[]
		zco=[]


		for j in range(titleline,len(inlines)):
			for string in elm:
				if inlines[j].find(string+"_")>-1 or inlines[j].find(string+" ")>-1 :
					if len(inlines[j].split()) > 3:
						atom.append(inlines[j])
		## The number of atoms 
		natom = len(atom)

		## Establish the molecular connectivity from the Gaussian input - This isn't important for attractive terms, however, it will dictate which repulsive terms are switched on
		conn=[]
		connectivity = 0

		for line in inlines:
			if "1 2 " in line:
				connectivity  = 1
				for j in range(natom+titleline+2,natom+titleline+2+natom):
					conn.append(inlines[j])

		bndidx=[]

		for j in range(0,natom):
			bndidx.append([0])
			for k in range(0,natom):
				bndidx[j].append(0)

		for j in range(0,natom):
			if connectivity == 1:
				for bonded in conn[j].split():
					if is_number(bonded) ==True:
						if int(bonded)-1!=j:
							bndidx[j][int(bonded)-1]=1
							bndidx[int(bonded)-1][j]=1

		#print atom
		for at in atom:
			atomtype.append((at.split()[0]).split("-")[0])
			xco.append(float(at.split()[1]))
			yco.append(float(at.split()[2]))
			zco.append(float(at.split()[3]))

		self.repulsive_vdw = 0.0
		self.attractive_vdw = 0.0
		
		
		molAatoms = getMollist(bndidx,0)
		mols = []
		for j in range(0,natom):
			mols.append(0)
			for atom in molAatoms:
				if atom == j: mols[j] = 1
		
		#print mols
		#print "o  Using the following C6 cooefficients"
		#print "   ID   Z    CN        C6"
		mxc=[0]
		for j in range(0,max_elem):
			mxc.append(0)
			for k in range(0,natom):
				if atomtype[k].find(elements[j])>-1:
					for l in range(0,maxc):
						if isinstance(c6ab[j][j][l][l], (list, tuple)):
							if c6ab[j][j][l][l][0]>0:
								#print "  ", atomtype[k],"  ", (j+1),"  ", c6ab[j][j][l][l][1],"  ", c6ab[j][j][l][l][0]
								mxc[j]=mxc[j]+1
					break

		## Coordination number based on covalent radii
		cn = ncoord(natom, rcov, atomtype, xco, yco, zco, max_elem, autoang, k1, k2)

		## C6 - Need to calculate these from fractional coordination
		#print "\n   #                XYZ [au]                   R0(AA) [Ang.]  CN          C6(AA)     C8(AA)   C10(AA) [au]"

		x=0
		for j in range(0,natom):
			C6jj = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,j)
			
			for k in range(0,natom):
				dum = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,k)
				x=x+dum
			#print j, k, x
			
			for k in range(0,max_elem):
				if atomtype[j].find(elements[k])>-1:z=k
			
			dum = 0.5*autoang*r[z][z]
			
			C8jj = 3.0*C6jj*math.pow(r2r4[z],2.0)
			C10jj=49.0/40.0 * math.pow(C8jj,2.0)/C6jj

		#	print "  ",(j+1), xco[j], yco[j], zco[j], atomtype[j], dum, cn[j], C6jj, C8jj, C10jj

		#print "\n   Molecular C6(AA) [au] =   ", x

		## Compute and output the individual components of the D3 energy correction ##
		#print "\n   Atoms  Types  C6            C8            E6              E8"
		for j in range(0,natom):
			
			## This could be used to 'switch off' dispersion between bonded or geminal atoms ##
			for k in range(j+1,natom):
				scalefactor=1.0
				vdw=1
				if bndidx[j][k]==1: vdw=0
				for l in range (0,natom):
					if bndidx[j][l] != 0 and bndidx[k][l]!=0 and j!=k and bndidx[j][k]==0: vdw=0
					for m in range (0,natom):
						if bndidx[j][l] != 0 and bndidx[l][m]!=0 and bndidx[k][m]!=0 and j!=m and k!=l and bndidx[j][m]==0: scalefactor=0.0
				
				
				if k>j and mols[j]!=mols[k]:
					#print "  ", (j+1), (k+1),"  ", atomtype[j], atomtype[k],"  ",
					
					## Pythagoras in 3D to work out distance ##
					xdist = xco[j]-xco[k]
					ydist = yco[j]-yco[k]
					zdist = zco[j]-zco[k]
					totdist = math.pow(xdist,2)+math.pow(ydist,2)+math.pow(zdist,2)
					totdist=math.sqrt(totdist)
					
					
					C6jk = getc6(maxc,max_elem,c6ab,mxc,atomtype,cn,j,k)
					
					## C8 parameters depend on C6 recursively
					for l in range(0,max_elem):
						if atomtype[j].find(elements[l])>-1:atomA=l
						if atomtype[k].find(elements[l])>-1:atomB=l
					
					C8jk = 3.0*C6jk*r2r4[atomA]*r2r4[atomB]
					C10jk=49.0/40.0 * math.pow(C8jk,2.0)/C6jk
					
					#print C6jk, C8jk,
					
					dist=totdist/autoang
					rr = r[atomA][atomB]/dist
					tmp1 = rs6*rr
					damp6 = 1/(1+6*math.pow(tmp1,alpha6))
					#print dist, r[atomA][atomB],rr, tmp1, damp6,
					tmp2 = rs8*rr
					damp8 = 1/(1+6*math.pow(tmp2,alpha8))
						
					## Repulsive coefficient
					R6jk = 4.5/math.pow(rs6,12)
					
					# Example for a repulsive potential dependent on R^-12
					# If two atoms are bonded then this term is zero
					if vdw == 1 and connectivity == 1:
						self.repulsive_vdw_term = R6jk*(1-damp6)/math.pow(1/tmp1,12)*math.pow(repfac,12)*scalefactor
					else: self.repulsive_vdw_term = 0.0
					
					#print R6jk, 1/math.pow(1/tmp1,12), math.pow(repfac,12)
					
					# Evaluation of the attractive term dependent on R^-6
					self.attractive_r6_term = -s6*C6jk*damp6/math.pow(dist,6)*autokcal*scalefactor
					
					# Evaluation of the attractive term dependent on R^-8
					self.attractive_r8_term = -s8*C8jk*damp8/math.pow(dist,8)*autokcal*scalefactor
					#self.attractive_r8_term = 0.0


					self.repulsive_vdw = self.repulsive_vdw + self.repulsive_vdw_term
					self.attractive_vdw = self.attractive_vdw + self.attractive_r6_term + self.attractive_r8_term
					#print self.attractive_r6_term, self.attractive_r8_term, self.repulsive_vdw_term
	
if __name__ == "__main__":
	
	# Takes arguments: (1) s6, (2) rs6, (3) s8, (4) input file(s)
	files = []
	if len(sys.argv) > 4:
		s6 = float(sys.argv[1])
		rs6 = float(sys.argv[2])
		s8 = float(sys.argv[3])
		for i in range(4,len(sys.argv)):
			files.append(sys.argv[i].split(".")[0])
	else:
		print "\nWrong number of arguments used. Correct format: eval_D3 s6 rs6 s8 file(s)\n"
		sys.exit()


	for file in files:
		fileD3 = calcD3(file, s6, rs6, s8)
		attractive_vdw = fileD3.attractive_vdw
		repulsive_vdw = fileD3.repulsive_vdw
		total_vdw = attractive_vdw + repulsive_vdw
		print "\no  Reading", file, "s6 =",s6, "rs6 = ", rs6, "s8 =",s8
		print "o  Attractive Part:", attractive_vdw, "kcal/mol"
		print "o  Repulsive Part:",repulsive_vdw, "kcal/mol"
		print "o  Total:", total_vdw,"kcal/mol \n"














