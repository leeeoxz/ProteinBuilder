import copy


class PDB_read:
	"""docstring for ClassName"""
	atom = "ATOM"
	alphaC = "CA"
	carb = ["C"]
	nitro = ["N"]
	endTag = "TER"
	modelTag = "MODEL"
	backboneAtom = ['N','CA','C','O']
	PDBname = None
	content = None
	atomList = None
	backbone = None
	coordinates = None
	alphaChain = None
	nChain = None
	cChain = None
	aaList = None

	def __init__(self, PDBname):
		self.PDBname = PDBname
		self.content = []
		with open (self.PDBname,"r") as f:
			for line in f:
				self.content.append(line)
		self.organizeAtom()
		self.setBackBone()
		self.setAlpha()
		self.setcChain()
		self.setnChain()
		self.setAA()

	def transforInt(self, atom): #transform coordinates and position into numbers
		if self.isNumber(atom[4]):
			for x,item in enumerate(atom):
				if x in range(5,10):
					atom[x] = float(atom[x]) #float for x,y,z coordinates
				elif x == 1:
					atom[x] = int(atom[x]) #float to atom number
				elif  x == 4:
					atom[x] = int(atom[x])
		else:
			for x,item in enumerate(atom):
				if x in range(6,11):
					atom[x] = float(atom[x]) #float for x,y,z coordinates
				elif x == 1:
					atom[x] = int(atom[x]) #float to atom number
				elif  x == 5:
					atom[x] = int(atom[x])
		return atom


	def organizeAtom(self):
		self.atomList = []
		for line in self.content:
			if self.atom in line[0:4]:
				line = line.replace("\n","").split()
				line = self.transforInt(line)
				self.atomList.append(line)
			elif self.endTag in line[0:3]:
				break

	def isNumber(self, value):
		try:
			int(value)
			return True

		except ValueError:
			return False

	def setBackBone(self):
		self.backbone = []
		for line in self.atomList:
			if line[2] in self.backboneAtom:
				self.backbone.append(line)

	def setAlpha(self):
		self.alphaChain = []
		for line in self.atomList:
			if self.alphaC in line[2]:
				self.alphaChain.append(line)

	def setCoord (self, atoms):
		self.coordinates =[]
		for lis in atoms:
			print atoms
			for x,item in enumerate(lis):
				print item,x
				if "." in item:
					atm = []
					atm.extend([float(lis[x]),float(lis[x+1]),float(lis[x+2])])
					self.coordinates.append(atm)
					break

	def setcChain(self):
		self.cChain = []
		for line in self.atomList:
			if line[2] in self.carb:
				self.cChain.append(line)

	def setAA(self):
		self.aaList = []
		for item in self.alphaChain:
			self.aaList.append(item[3])

	def setnChain(self):
		self.nChain = []
		for line in self.atomList:
			if line[2] in self.nitro:
				self.nChain.append(line)

	def getAtom(self):
		return copy.deepcopy(self.atomList)

	def getBackBone(self):
		return copy.deepcopy(self.backbone)

	def getCoord (self,atom):
		self.setCoord(atom)
		return copy.deepcopy(self.coordinates)

	def getAlpha(self):
		return copy.deepcopy(self.alphaChain)

	def getN(self):
		return copy.deepcopy(self.nChain)

	def getC(self):
		return copy.deepcopy(self.cChain)

	def getAA(self):
		return copy.deepcopy(self.aaList)
