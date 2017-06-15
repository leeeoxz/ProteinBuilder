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
	transRot = None #Translational and rotational vector

	def __init__(self,file1,file2, parameter):
		self.pdb1 = PDB_read(file1) #Opens the first file as reference
		self.pdb2 = PDB_read(file2) #Opens the second file to be the one compared
		self.sm = SimAnne(6,[-100.0,-100.0,-100.0,0.0,0.0,0.0],[100.0,100.0,100.0,360.0,360.0,360.0]) #Dimension, translocation ranges, rotation ranges
		self.parameter = parameter #Defines wich structure part will be used (AlphaChain,AllAtoms,Backbone)
		self.transRot = self.sm.create_solution()
		


	
	def Rotate(self, transRot,atomList):
		#rotate around x axys
		xRotation = np.matrix([[1.0,0.0,0.0],[0,math.cos(math.radians(transRot[3])),(-math.sin(math.radians(transRot[3])))],[0,math.sin(math.radians(transRot[3])),math.cos(math.radians(transRot[3]))]])
		#rotate around y axys
		yRotation = np.matrix([[math.cos(math.radians(transRot[4])),0.0,math.sin(math.radians(transRot[4]))],[0.0,1.0,0.0],[(-math.sin(math.radians(transRot[4]))),0.0,math.cos(math.radians(transRot[4]))]])
		#rotate around z axys
		zRotation = np.matrix([[math.cos(math.radians(transRot[5])),(-math.sin(math.radians(transRot[5]))),0.0],[math.sin(math.radians(transRot[5])),math.cos(math.radians(transRot[5])),0.0],[0.0,0.0,1.0]])
		Rotation = xRotation*zRotation*yRotation
		newAtom = atomList * Rotation


	def Translate (self,transRot, atomList):
		#Transform the translational factor to a matrix, due to math lib functions
		transMatrix=np.matrix([transRot[0:3]]*len(atomList))
		#translocation:
		newAtom = atomList + transMatrix
		return newAtom

	def setRMSD (self,cordRef, cordPDB, TransRot):
		aligned = Translate(TransRot,cordPDB)
		aligned = Rotate(aligned,cordPDB)
		rm = aligned - cordRef
		rm = np.power(rm,2)
		rm = np.sum(rm)
		rm = math.sqrt(float(rm))
		rm = rm/float(len(cordRef))
		return rm

	def rotine(self):
		while self.sm.getTemperature() > 1: #minimum temperature
			iterations = 2000
			while iterations > 0:
				newTransRot = map(float,self.sm.create_solution())
				solution1 = self.setRMSD(self.pdb1,self.pdb2,newTransRot)
				if solution1 - self.bestSolution <=0:
					self.bestSolution = solution1
					self.transRot = newTransRot
				elif self.sm.getProbability(solution1,self.bestSolution,self.sm.getTemperature()) > random.uniform(0.0, 1.0):
					self.bestSolution = solution1
					self.transRot = newTransRot
				iterations -= 1
			sm.updateTemperature()