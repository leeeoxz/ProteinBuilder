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

atomsN =["N","H1","H2","H3"]
atomsCO = ["C","O"]

def rotate (angle, pos): #ROTATES DE COORDINATES
	xRotation = np.matrix([[1.0,0.0,0.0],[0,math.cos(math.radians(angle)),(-math.sin(math.radians(angle)))],[0,math.sin(math.radians(angle)),math.cos(math.radians(angle))]])
	#rotate around y axys
	yRotation = np.matrix([[math.cos(math.radians(angle)),0.0,math.sin(math.radians(angle))],[0.0,1.0,0.0],[(-math.sin(math.radians(angle))),0.0,math.cos(math.radians(angle))]])
	#rotate around z axys
	zRotation = np.matrix([[math.cos(math.radians(angle)),(-math.sin(math.radians(angle))),0.0],[math.sin(math.radians(angle)),math.cos(math.radians(angle)),0.0],[0.0,0.0,1.0]])
	Rotation = xRotation*zRotation*yRotation
	newAtom = pos * Rotation
	return newAtom

def psi(structure,angle):

	for x,aa in enumerate(structure):
		if x == 0:
			for atom in aa:
				print atom
				if atom[2] in atomsCO:
					print atom
					atom[5:8] = rotate(angle, atom[5:8])
					print atom

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

aaNuber = len(coordAA1) #GETS THE NUMBER OF AMINO ACIDS
tam = (aaNuber*2)-2 #SET THE NUMBERS OF ANGLES 


sm = SimAnne(tam,[0.0]*tam,[360.0]*tam) # STARTS THE METAHEURISTICS

angles = sm.create_solution()
for x,aa in enumerate(coordAA1):
	if x == 0:
		psi(coordAA1,angles[0])