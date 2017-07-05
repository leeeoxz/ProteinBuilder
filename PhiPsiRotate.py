# Leonardo Santos #
# June 2017 #

# Tries to find the best dihedral angles to a structure #

from PDB_read import PDB_read
from SA import SimAnne
import rmsd
import math
import numpy as np
import copy
import random
from matplotlib import pyplot as plt


atomsPHI =["N","H1","H2","H3","1H","2H","3H","CA"]
atomsPSI = ["C"]

def calc_exact_rmsd(ref_atoms, mob_atoms):
    ra = np.array(copy.deepcopy(ref_atoms))
    ma = np.array(copy.deepcopy(mob_atoms))
    ra -= rmsd.centroid(ra)
    ma -= rmsd.centroid(ma)
    return rmsd.kabsch_rmsd(ra, ma)

def rotate(angle,atm,carb,nitro):
	ang = angle
	v = np.array(atm) - np.array(carb)
	k = np.array(nitro) - np.array(carb)
	k = k/np.linalg.norm(k)
	v = v * np.cos(ang) + (np.cross(k,v))*np.sin(ang) + k*(np.dot(k,v)) *(1.0-np.cos(ang))
	v = v + np.array(carb)
	return list(v)

def psi(structure,angle):
	newStruct = []
	for atom in structure[0]:
		if atom[2] == "CA":
			pCA = atom[5:8]
		elif atom[2] == "C":
			pC = atom[5:8]
	for x,aa in enumerate(structure):
		amino = []
		if x == 0:
			for atom in aa:
				if atom[2]== "O":
					atom[5:8] = rotate(angle, atom[5:8],pCA,pC)
					amino.append(atom)
				else:
					amino.append(atom)
		else:
			for atom in aa:
				atom[5:8] = rotate(angle, atom[5:8],pCA,pC)
				amino.append(atom)
		newStruct.append(amino)
	return newStruct

def phi(structure,angle):
	newStruct = []
	for atom in structure[0]:
		if atom[2] == "CA":
			pCA = atom[5:8]
		elif atom[2] == "N":
			pN = atom[5:8]
	for x,aa in enumerate(structure):
		amino = []
		if x == 0:
			for atom in aa:
				if atom[2] not in atomsPHI:
					atom[5:8] = rotate(angle, atom[5:8],pN,pCA)
					amino.append(atom)
				else:
					amino.append(atom)
		else:
			for atom in aa:
				atom[5:8] = rotate(angle, atom[5:8],pN,pCA)
				amino.append(atom)
		newStruct.append(amino)
	return newStruct

def diff(target, ang):
	a = ang - target
	if a > math.pi:
		a -= math.pi *2.
	if a < -math.pi:
		a += math.pi *2.
	return a

def createAngles():
	anglesPHI = sm.create_solution()
	anglesPSI = sm.create_solution()
	for x,angle in enumerate(anglesPHI):
		anglesPHI[x] = diff(math.radians(180.),math.radians(angle))
	for x,angle in enumerate(anglesPSI):
		anglesPSI[x] = diff(math.radians(180.),math.radians(angle))
	angles = zip(anglesPHI,anglesPSI)
	return angles

def withoutAA(atm):
	v =[]
	for aa in atm:
		for atom in aa:
			v.append(atom)
	return v

def isNumber(value):
	try:
		int(value)
		return True

	except ValueError:
		return False

#angles =  [[diff(math.radians(360.0),math.radians(360.0)),diff(math.radians(180.),math.radians(176.63))],[diff(math.radians(180.),math.radians(148.48)),diff(math.radians(180.),math.radians(-21.96))],[diff(math.radians(180.),math.radians(114.02)),diff(math.radians(180.),math.radians(29.89))],[diff(math.radians(180.),math.radians(-88.0)),diff(math.radians(180.),math.radians(-38.16))],[diff(math.radians(180.),math.radians(-74.24)),diff(math.radians(180.),math.radians(360.0))]]

def AAsetter(atoms,tam):
	atomsAA = []
	if isNumber(atoms[0][4]):
		index = 4
	else:
		index = 5
	for i in xrange(tam):
		aa = []
		for atom in atoms:
			if atom[index] == i+1:
				aa.append(atom)
		atomsAA.append(aa)
	return atomsAA


def rotatephipsi(atoms,angles):
	for x,aa in enumerate(atoms):
		if x == 0:
			atoms = psi(atoms,angles[x][1])
		elif x == (tam-1):
			atoms = atoms[:x]+(phi(atoms[x:],angles[x][0]))
		else:
			atoms = atoms[:x]+(phi(atoms[x:],angles[x][0]))
			atoms = atoms[:x]+(psi(atoms[x:],angles[x][1]))
	return atoms

file1=PDB_read("YGGFM.pdb") 

atoms1 = file1.getAtomsAA() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST
atoms2 = PDB_read("1plx.pdb").getAtomsAA() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST

coordAA1 =[]
coordAA2 =[]

for i in xrange(len(atoms1)):
	coordAA1.append(file1.getCoord(atoms1[i]))
for i in xrange(len(atoms2)):
	coordAA2.append(PDB_read("1plx.pdb").getCoord(atoms2[i]))


allatoms2 = []
allatoms1 = []
for item in coordAA2:
	for i in item:
		allatoms2.append(i)
for item in coordAA1:
	for i in item:
		allatoms1.append(i)


s0 = calc_exact_rmsd(allatoms2,allatoms1)

tam = len(coordAA1) #GETS THE NUMBER OF AMINO ACIDS

sm = SimAnne(tam,[0.0]*tam,[360.0]*tam) # STARTS THE METAHEURISTICS

ang = createAngles() #INITIAL SOLUTION

atoms1 = rotatephipsi(atoms1,ang) # ROTATE PHI AND PSI WITH THE FIRST SOLUTION
#atoms1 = withoutAA(atoms1)
#atoms1 = file1.getCoord(atoms1)

best = s0
bestst = atoms1
bestn = 0
n = 0 
m = []
p = []
while sm.getTemperature() > 1:
	iterations = 100
	while iterations > 0:
		ang = createAngles()
		newTransRot = rotatephipsi(atoms1,ang)
		newTransRot = withoutAA(newTransRot)
		allatoms1 = file1.getCoord(newTransRot)
		s1 = calc_exact_rmsd(allatoms2,allatoms1)
		if s1<s0:
			s0 = s1
			atoms1 = AAsetter(newTransRot,tam)
			if s0 < best:
				best = s0
				bestst = atoms1
				bestn = sm.getTemperature()
		elif sm.getProbability(s1,s0,sm.getTemperature()) > random.uniform(0.0, 1.0):
			s0 = s1
			atoms1 = AAsetter(newTransRot,tam)
		iterations -= 1
		n += 1
		m.append(n)
		p.append(s0)	
	sm.updateTemperature()


file = open("1PLX-F.pdb","w") #Opens the output file

for aa in bestst:
	for atom in aa:
		line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(atom[0], atom[1], atom[2], " ", atom[3], " ", atom[4], " ", atom[5], atom[6], atom[7], atom[8], atom[9], " ", " ")
		file.write(line)
		file.write("\n")
file.write("TER")
file.close()

plt.plot(m,p)
plt.xlim(len(m), 0)  # decreasing time

plt.xlabel('RMSD')
plt.ylabel('Temperature')
plt.title('YGGFM; Best RMSD: '+str(best))
plt.savefig('destination_path.eps', format='eps', dpi=1000)
plt.show()