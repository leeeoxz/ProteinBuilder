from PDB_read import PDB_read
from SA import SimAnne
import math
import numpy as np
import copy
import random


class rmsd():
	pdb1 = None #Reference structure; .pdb file 
	pdb2 = None #Structure to be compared; .pdb file
	sm = None #Simulated Annealing class starter
	bestSolution = None #Best RMSD

	def __init__(self,file1,file2, parameter):
		self.pdb1 = PDB_read(file1) #Opens the first file as reference
		self.pdb2 = PDB_read(file2) #Opens the second file to be the one compared
		self.sm = SimAnne(6,[-100.0,-100.0,-100.0,0.0,0.0,0.0],[100.0,100.0,100.0,360.0,360.0,360.0]) #Dimension, translocation ranges, rotation ranges
		self.parameter = parameter #Defines wich structure part will be used (AlphaChain,AllAtoms,Backbone)
		self.bestSolution = self.calcRMSD() 




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