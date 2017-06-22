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

def rotate (angle, pos): #ROTATES DE COORDINATES
	xRotation = np.matrix([[1.0,0.0,0.0],[0,math.cos(math.radians(angle)),(-math.sin(math.radians(angle)))],[0,math.sin(math.radians(angle)),math.cos(math.radians(angle))]])
	#rotate around y axys
	yRotation = np.matrix([[math.cos(math.radians(angle)),0.0,math.sin(math.radians(angle))],[0.0,1.0,0.0],[(-math.sin(math.radians(angle))),0.0,math.cos(math.radians(angle))]])
	#rotate around z axys
	zRotation = np.matrix([[math.cos(math.radians(angle)),(-math.sin(math.radians(angle))),0.0],[math.sin(math.radians(angle)),math.cos(math.radians(angle)),0.0],[0.0,0.0,1.0]])
	Rotation = xRotation*zRotation*yRotation
	newAtom = pos * Rotation
	a = np.array(newAtom)[0]
	return a

def psi(structure,angle):
	newStruct = []
	for x,aa in enumerate(structure):
		amino = []
		if x == 0:
			for atom in aa:
				if atom[2] in atomsPSI:
					atom[5:8] = rotate(angle, atom[5:8])
					amino.append(atom)
				else:
					amino.append(atom)
		else:
			for atom in aa:
				atom[5:8] = rotate(angle, atom[5:8])
				amino.append(atom)
		newStruct.append(amino)
	return newStruct

def phi(structure,angle):
	newStruct = []
	for x,aa in enumerate(structure):
		amino = []
		if x == 0:
			for atom in aa:
				if atom[2] not in atomsPHI:
					atom[5:8] = rotate(angle, atom[5:8])
					amino.append(atom)
				else:
					amino.append(atom)
		else:
			for atom in aa:
				atom[5:8] = rotate(angle, atom[5:8])
				amino.append(atom)
		newStruct.append(amino)
	return newStruct


# TEM QUE MANDAR TODO O ATOMO, NAO SO AS COORDENADAS FUCK #

file1=PDB_read("YGGFM.pdb") 

atoms1 = file1.getAtomsAA() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST
atoms2 = PDB_read("1plx.pdb").getAtomsAA() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST
coordAA1 =[]
coordAA2 =[]
for i in xrange(len(atoms1)):
	coordAA1.append(file1.getCoord(atoms1[i]))
for i in xrange(len(atoms2)):
	coordAA2.append(PDB_read("1plx.pdb").getCoord(atoms2[i]))

tam = len(coordAA1) #GETS THE NUMBER OF AMINO ACIDS


sm = SimAnne(tam,[0.0]*tam,[360.0]*tam) # STARTS THE METAHEURISTICS

anglesPHI = sm.create_solution()
anglesPSI = sm.create_solution()

angles = zip(anglesPHI,anglesPSI)



for x,aa in enumerate(atoms1):
	if x == 0:
		atoms1 = psi(atoms1,angles[x][0])
	elif x == (tam-1):
		atoms1 = atoms1[:x]+(phi(atoms1[x:],angles[x][0]))
	else:
		atoms1 = atoms1[:x]+(phi(atoms1[x:],angles[x][0]))
		atoms1 = atoms1[:x]+(psi(atoms1[x:],angles[x][1]))

print atoms1
file = open("1PLX-F.pdb","w") #Opens the output file

for aa in atoms1:
	for atom in aa:
		line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(atom[0], atom[1], atom[2], " ", atom[3], " ", atom[4], " ", atom[5], atom[6], atom[7], atom[8], atom[9], " ", " ")
		file.write(line)
		file.write("\n")
file.write("TER")
file.close()