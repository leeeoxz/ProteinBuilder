# Leonardo Santos #
# June 2017 #

# Tries to find the best dihedral angles to a structure #

from PDB_read import PDB_read
from SA import SimAnne
from RMSD import rmsd
import math
import numpy as np
import copy
import random

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
					print atom[2]
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


# TEM QUE MANDAR TODO O ATOMO, NAO SO AS COORDENADAS FUCK #

file1=PDB_read("NLYIQWLKDGGPSSGRPPPS.pdb") 

atoms1 = file1.getAtomsAA() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST
#atoms2 = PDB_read("1plx.pdb").getAtomsAA() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST
coordAA1 =[]
coordAA2 =[]
for i in xrange(len(atoms1)):
	coordAA1.append(file1.getCoord(atoms1[i]))
#for i in xrange(len(atoms2)):
#	coordAA2.append(PDB_read("1plx.pdb").getCoord(atoms2[i]))

tam = len(coordAA1) #GETS THE NUMBER OF AMINO ACIDS


sm = SimAnne(tam,[0.0]*tam,[360.0]*tam) # STARTS THE METAHEURISTICS

anglesPHI = sm.create_solution()
anglesPSI = sm.create_solution()

#angles = zip(anglesPHI,anglesPSI)
angles =  [[diff(math.radians(360.0),math.radians(360.0)),diff(math.radians(180.),math.radians(0.00))],[diff(math.radians(180.),math.radians(-172.67)),diff(math.radians(180.),math.radians(-62.14))],[diff(math.radians(180.),math.radians(-151.94)),diff(math.radians(180.),math.radians(180.00))],[diff(math.radians(180.),math.radians(132.87)),diff(math.radians(180.),math.radians(0.12))],[diff(math.radians(180.),math.radians(68.62)),diff(math.radians(180.),math.radians(0.01))],
[diff(math.radians(180.),math.radians(-42.86)),diff(math.radians(180.),math.radians(-179.99))],[diff(math.radians(180.),math.radians(-8.77)),diff(math.radians(180.),math.radians(-62.07))],[diff(math.radians(180.),math.radians(74.85)),diff(math.radians(180.),math.radians(-179.98))],[diff(math.radians(180.),math.radians(122.09)),diff(math.radians(180.),math.radians(-178.84))],[diff(math.radians(180.),math.radians(88.06)),diff(math.radians(180.),math.radians(179.97))],
[diff(math.radians(180.),math.radians(80.46)),diff(math.radians(180.),math.radians(-179.94))],[diff(math.radians(180.),math.radians(107.72)),diff(math.radians(180.),math.radians(179.94))],[diff(math.radians(180.),math.radians(116.51)),diff(math.radians(180.),math.radians(-179.99))],[diff(math.radians(180.),math.radians(99.99)),diff(math.radians(180.),math.radians(179.99))],[diff(math.radians(180.),math.radians(159.44)),diff(math.radians(180.),math.radians(179.99))],
[diff(math.radians(180.),math.radians(66.67)),diff(math.radians(180.),math.radians(178.77))],[diff(math.radians(180.),math.radians(131.47)),diff(math.radians(180.),math.radians(179.98))],[diff(math.radians(180.),math.radians(-127.43)),diff(math.radians(180.),math.radians(179.97))],[diff(math.radians(180.),math.radians(-170.73)),diff(math.radians(180.),math.radians(-179.99))],[diff(math.radians(180.),math.radians(-7.24)),diff(math.radians(360.),math.radians(360.))],]
print angles

for x,aa in enumerate(atoms1):
	if x == 0:
		atoms1 = psi(atoms1,angles[x][1])
	elif x == (tam-1):
		atoms1 = atoms1[:x]+(phi(atoms1[x:],angles[x][0]))
	else:
		atoms1 = atoms1[:x]+(phi(atoms1[x:],angles[x][0]))
		atoms1 = atoms1[:x]+(psi(atoms1[x:],angles[x][1]))


file = open("1L2Y-P.pdb","w") #Opens the output file

for aa in atoms1:
	for atom in aa:
		line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(atom[0], atom[1], atom[2], " ", atom[3], " ", atom[4], " ", atom[5], atom[6], atom[7], atom[8], atom[9], " ", " ")
		t = line[11:17]
		a = atom[-1]
		print a
		line = line.replace(t,a)
		file.write(line)
		file.write("\n")
file.write("TER")
file.close()
