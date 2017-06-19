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

atoms1 = PDB_read("1plx.pdb").getAtom() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST
atoms2 = PDB_read("YGGFM.pdb").getAtom() # GET THE INFORMATIONS OF ALL ATOMS AS A LIST
coord1 = PDB_read("1plx.pdb").getCoord(atoms1) # GET THE COORDENATES OF ALL ATOMS
coord2 = PDB_read("YGGFM.pdb").getCoord(atoms2) # GET THE COORDENATES OF ALL ATOMS

rms = rmsd(coord1,coord2,"backbone") 

tam = (len(PDB_read("YGGFM.pdb").getAA())*2)-2 #GETS THE NUMBER OF AMINO ACIDS

sm = SimAnne(tam,[0.0]*tam,[360.0]*tam) # STARTS THE METAHEURISTICS


def rotate (angle, pos): #ROTATES DE COORDINATES
	xRotation = np.matrix([[1.0,0.0,0.0],[0,math.cos(math.radians(angle)),(-math.sin(math.radians(angle)))],[0,math.sin(math.radians(angle)),math.cos(math.radians(angle))]])
	#rotate around y axys
	yRotation = np.matrix([[math.cos(math.radians(angle)),0.0,math.sin(math.radians(angle))],[0.0,1.0,0.0],[(-math.sin(math.radians(angle))),0.0,math.cos(math.radians(angle))]])
	#rotate around z axys
	zRotation = np.matrix([[math.cos(math.radians(angle)),(-math.sin(math.radians(angle))),0.0],[math.sin(math.radians(angle)),math.cos(math.radians(angle)),0.0],[0.0,0.0,1.0]])
	Rotation = xRotation*zRotation*yRotation
	newAtom = pos * Rotation
	return newAtom

def setStart(value):
	try:
		int(value[4])
		return 5

	except ValueError:
		return 6

start = setStart(atoms2[0])

for x,atom in enumerate(atoms2):
	if x == 0 : 
		pass		