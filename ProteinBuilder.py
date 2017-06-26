from AADictBuilder import BuildDictionary as BD
import copy
from Rotational import getAngle
import numpy as np
import math

class ProteinBuilder():
	aminoacidDict = None
	carbTAG = ["OC","HOC","HO","HC"]
	hTAG = ["1H","H"]
	seq = None
	proteinStruct = None #array of amino acids, which are arrays of atoms
	seqSize = None #int value representing the len of the sequence

	def __init__(self,sequence):
		self.seq = sequence
		self.aminoacidDict = BD().getDic() #creates a dictionary with all 20 aminoacids and its atom locations
		self.proteinStruct = []
		self.setStructure()
		self.seqSize = len(self.proteinStruct)
		self.setAngle()


	def removeOH(self,last):
		index = []
		for ii,atom in enumerate(last):
			if atom[2] in self.carbTAG: #if its O or H from OH it saves:
				index.append(ii) #the index, so we can delete it, and
		for i in sorted(index,reverse=True): #delete the atoms
			del last[i]
		return last

	def removeH(self, new):
		for x,atom in enumerate(new):
			if atom[2] in self.hTAG:
				index = x
		del new[index]
		return new

	def getDist (self,last, new):  #last and new are amino acids, that is two lists
		for ii,atom in enumerate(last):
			if atom[2] in self.carbTAG:
				x = atom[5] - new[0][5] #distance to move on x axis
				y = atom[6]	- new[0][6]	#distance to move on y axis
				z = atom[7]	- new[0][7] #distance to move on z axis
				n_aa = atom[4] + 1 #get the number of amino acid
				return n_aa,x,y,z


	def addEnd (self,structure, aa):
		dist= self.getDist(structure[-1],aa) #gets a list (number of next atom, number of next amino acid, x distance,y distance,z distance)
		structure[-1] = self.removeOH(structure[-1])
		aa = self.removeH(aa)
		cont = 1
		for atom in aa:
			atom[4] = dist[0] #set the new position of the aa in sequence
			atom[5] = round(atom[5] + dist[1],3) #translocate the aminoacid on x axis
			atom[6] = round(atom[6] + dist[2],3) #translocate the aminoacid on y axis
			atom[7] = round(atom[7] + dist[3],3) #translocate the aminoacid on z axis
		structure.append(aa)
			
	def setStructure(self):
		for x,aa in enumerate(self.seq):
			if not self.proteinStruct: #if proteinStructu is empty than its the first amino acid
				self.proteinStruct.append(copy.deepcopy(self.aminoacidDict[aa]))
			else:
				self.addEnd(self.proteinStruct,copy.deepcopy(self.aminoacidDict[aa]))
		x = 1
		for aa in self.proteinStruct: # set the atom's numbers
			for atom in aa:
				atom[1] = x
				x = x +1

	def rotate(self):
		

	def calcAngle3p(self,ap,cp,pp):
		ap_mt = np.matrix(ap)
		cp_mt = np.matrix(cp)
		pp_mt = np.matrix(pp)
		v1 = ap_mt - cp_mt
		v2 = cp_mt - pp_mt
		v1n= v1/np.linalg.norm(v1)
		v2n= v2/np.linalg.norm(v2)
		res = np.sum(np.multiply(v1n,v2n))
		return math.degrees(float(np.arccos(res)))

	def setAngle(self):
		for x,aa in enumerate(self.proteinStruct):
			if (x!=(self.seqSize-1)): #first amino acid
				for atom in aa:
					if atom[2] == 'C': #only calculates from C terminal as the anterior position
						ap = [atom[5],atom[6],atom[7]] #ap = anterior position, [float(Px),float(Py),float(Pz)]
						for atom in self.proteinStruct[x+1]: #gets the coordinates from the next amino acid
							if atom[2] == "N": #in this case, the Nitrogen is going to be the central atom
								cp = [atom[5],atom[6],atom[7]] #cp = central position, [float(Px),float(Py),float(Pz)]
							elif atom[2] == "2H":  #in this case, the Hydrogen bonded to Nitrogen is going to be the posterial atom
							    ppH = [atom[5],atom[6],atom[7]] #pp = posteral position, [float(Px),float(Py),float(Pz)]
							elif atom[2] == "CA": #in this case, the Alpha Carbon bonded to Nitrogen is going to be the posterial atom
								ppCA = [atom[5],atom[6],atom[7]] #pp = posteral position, [float(Px),float(Py),float(Pz)]
			elif (x==(self.seqSize-1)):
				break
			rotate = 120 - self.calcAngle3p(ap,cp,ppH)
			print rotate



	def getStructure (self):
		return copy.deepcopy(self.proteinStruct)


ProteinBuilder("YGGMF")