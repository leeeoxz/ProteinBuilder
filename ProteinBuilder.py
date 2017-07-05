from AADictBuilder import BuildDictionary as BD
import copy
from Rotational import getAngle as getAngle
import numpy as np
import math

class ProteinBuilder():
	aminoacidDict = None
	carbTAG = ["OC","HOC","HO","HC"]
	hTAG = ["2H","H"]
	angleH = ["1H"]
	alphaC = ["CA"]
	nitro = ["N"]
	cabTAG = ["C"]
	seq = None
	proteinStruct = None #array of amino acids, which are arrays of atoms
	seqSize = None #int value representing the len of the sequence

	def __init__(self,sequence):
		self.seq = sequence
		self.aminoacidDict = BD().getDic() #creates a dictionary with all 20 aminoacids and its atom locations
		self.proteinStruct = []
		self.setStructure()
		self.seqSize = len(self.proteinStruct)
		#self.setAngle()


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
		for x,atom in enumerate(structure[-1]):
			if atom[2] in self.cabTAG:
				paC = atom[5:8]
			elif atom[2] in self.alphaC:
				paCA = atom[5:8]
		structure[-1] = self.removeOH(structure[-1])
		aa = self.removeH(aa)
		cont = 1
		for x,atom in enumerate(aa):
			atom[4] = dist[0] #set the new position of the aa in sequence
			atom[5] = round(atom[5] + dist[1],3) #translocate the aminoacid on x axis
			atom[6] = round(atom[6] + dist[2],3) #translocate the aminoacid on y axis
			atom[7] = round(atom[7] + dist[3],3) #translocate the aminoacid on z axis
			if atom[2] in self.angleH: #get the coordinates of H to rotate
				H = x
				ppH = atom[5:8]
			elif atom[2] in self.alphaC: #get the coordinates of H to rotate
				CA = x
				ppCA = atom[5:8]
			elif atom[2] in self.nitro: #get the coordinates of H to rotate
				ppN = atom[5:8]
				N = x

		if aa[0][3] != "PRO":
			angle = self.calcAngle3p(paC,ppN,ppH)
			angH =  self.diff(120., angle[0])
			if angH > angle[0]:
				angH = self.diff(120., -angle[0])
	 
			ortho_aux = np.ndarray.tolist(angle[1])[0]
			ortho = []
			for item in ortho_aux:
				ortho.append(round(item))
			aa[H][5:8] = self.rotate(ortho,angH,ppH,paC,ppN,ppH)

			
			angle = self.calcAngle3p(paC,ppN,ppCA)
			angCA = self.diff(120.,angle[0])
			if angCA > angle[0]:
				angCA = self.diff(120., -angle[0])

			ortho_aux = np.ndarray.tolist(angle[1])[0]
			ortho = []
			for item in ortho_aux:
				ortho.append(round(item))
			
			for x,atom in enumerate(aa[1:]):
				if aa[x+1][2] not in self.angleH:
					aa[x+1][5:8] = self.rotate(ortho,angCA,aa[x+1][5:8],paC,ppN,ppCA)
		
		omega = getAngle.getAngle(getAngle(paCA,paC,ppN,ppCA))
		torotate = self.diff(180.,omega)
		for x,atom in enumerate(aa):
			aa[x][5:8] =self.rotateOmega(torotate,aa[x][5:8],paC,ppN)
		
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

	def diff(self,target, ang):
		a = math.radians(target) - math.radians(ang)
		if a > math.pi:
			a -= math.pi *2.
		if a < -math.pi:
			a += math.pi *2.
		return a


	def rotate(self,ortho,ang,atm,at,ct,ps):
		cos = np.cos(ang)
		sen = np.sin(ang)
		tan = 1.0-cos
		pos = np.array(atm) - np.array(ct)
		rot = np.matrix([[cos+ortho[0]*ortho[0]*tan, ortho[0] * ortho[1] * tan - ortho[2] * sen, ortho[0] * ortho[2] * tan + ortho[1] * sen],[ortho[0] * ortho[1] * tan + ortho[2] * sen, cos + ortho[1] * ortho[1] * tan, ortho[1] * ortho[2] * tan - ortho[0] * sen], [ortho[2] * ortho[0] * tan - ortho[1] * sen, ortho[2] * ortho[1] * tan + ortho[0] * sen, cos + ortho[2] * ortho[2] * tan]])
		new_pos = np.matrix.tolist(np.matrix(pos) * rot.transpose())[0]
		new_pos = list(np.array(new_pos) + np.array(ct))
		return new_pos

	def calcAngle3p(self,ap,cp,pp):
		ap_mt = np.matrix(ap)
		cp_mt = np.matrix(cp)
		pp_mt = np.matrix(pp)
		v1 = ap_mt - cp_mt
		v2 = pp_mt - cp_mt
		v1n= v1/np.linalg.norm(v1)
		v2n= v2/np.linalg.norm(v2)
		res = np.sum(np.multiply(v1n,v2n))
		return math.degrees(float(np.arccos(res))),np.cross(v1,v2)/np.linalg.norm(np.cross(v1,v2))

	def rotateOmega(self,angle,atm,carb,nitro):
		ang = angle
		v = np.array(atm) - np.array(carb)
		k = np.array(nitro) - np.array(carb)
		k = k/np.linalg.norm(k)
		v = v * np.cos(ang) + (np.cross(k,v))*np.sin(ang) + k*(np.dot(k,v)) *(1.0-np.cos(ang))
		v = v + np.array(carb)
		return list(v)

	def getStructure (self):
		return copy.deepcopy(self.proteinStruct)