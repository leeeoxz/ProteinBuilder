import copy
import numpy as np
import math


class getAngle():
	b1 = None #Distance vector between atoms 1,2
	b2 = None #Distance vector between atoms 2,3
	b3 = None #Distance vector between atoms 3,4
	atoms = None #Gets vectors with atom positions, if calculates phi gets vector with [C1,N2,CA2,C2], if calculates psi gets vector with [N2,CA2,C2,N3]
	n1 = None 
	n2 = None
	angle = None #desired angle
	b2V = None
	m1 = None
	x = None
	y = None

	def __init__(self,a1,a2,a3,a4):
		self.atoms = [np.array(a1),np.array(a2),np.array(a3),np.array(a4)]
		self.b1 = self.atoms[0] - self.atoms[1]
		self.b2 = self.atoms[1] - self.atoms[2]
		self.b3 = self.atoms[2] - self.atoms[3]
		self.n1 = self.calcN(self.b1,self.b2)
		self.n2 = self.calcN(self.b2,self.b3)
		self.calcB2V()
		self.calcM()
		self.x = np.inner(self.n1,self.n2)
		self.y = np.inner(self.m1,self.n2)
		self.angle = round( math.degrees( math.atan2( self.y, self.x ) ), 2 )

	def calcN(self,bx,by):
		i = np.cross(bx,by)
		j = np.linalg.norm(i)
		return (i/j)

	def calcM(self):
		self.m1 = np.cross(self.n1,self.b2V)

	def calcB2V(self):
		self.b2V = self.b2/(np.linalg.norm(self.b2))

	def getAngle(self):
		return copy.deepcopy(self.angle)

