from PDB_read import PDB_read
from SA import SimAnne
import math
import numpy as np
import copy
import random

def aligner (transRot, atomList):
	xRotation = np.matrix([[1.0,0.0,0.0],[0,math.cos(math.radians(transRot[3])),(-math.sin(math.radians(transRot[3])))],[0,math.sin(math.radians(transRot[3])),math.cos(math.radians(transRot[3]))]])
	yRotation = np.matrix([[math.cos(math.radians(transRot[4])),0.0,math.sin(math.radians(transRot[4]))],[0.0,1.0,0.0],[(-math.sin(math.radians(transRot[4]))),0.0,math.cos(math.radians(transRot[4]))]])
	zRotation = np.matrix([[math.cos(math.radians(transRot[5])),(-math.sin(math.radians(transRot[5]))),0.0],[math.sin(math.radians(transRot[5])),math.cos(math.radians(transRot[5])),0.0],[0.0,0.0,1.0]])
	Rotation = xRotation*zRotation*yRotation
	#Transform the translational factor to a matrix, due to math lib functions
	transMatrix=np.matrix([transRot[0:3]]*len(atomList))
	#translocation:
	newAtom = atomList + transMatrix
	#rotation:
	newAtom = newAtom * Rotation
	return newAtom



def rmsd (cordRef, cordPDB, TransRot):
	aligned = aligner(TransRot,cordPDB)
	rm = aligned -cordRef
	rm = np.power(rm,2)
	rm = np.sum(rm)
	rm = math.sqrt(float(rm))
	rm = rm/float(len(cordRef))
	return rm

#reference = raw_input("Reference input: ")
#pdb = raw_input("PDB file: ")

ref = PDB_read("reference.pdb").getCoord()
pdbComp = PDB_read("1ACW-02.pdb").getCoord()


simu_ann = SimAnne(6,[-100.0,-100.0,-100.0,0.0,0.0,0.0],[100.0,100.0,100.0,360.0,360.0,360.0])

simu_ann.solution = map(float,simu_ann.solution)

s0 = rmsd(ref,pdbComp,simu_ann.solution)

while simu_ann.getTemperature() > 1:
	iterations = 2000
	while iterations > 0:
		newTransRot = map(float,simu_ann.create_solution())
		s1 = rmsd(ref,pdbComp,newTransRot)
		if s1-s0 <= 0:
			s0 = s1
			v = newTransRot
		elif simu_ann.getProbability(s1,s0,simu_ann.getTemperature()) > random.uniform(0.0, 1.0):
			s0 = s1
			v = newTransRot
		iterations -= 1
	simu_ann.updateTemperature()


print s0, v