from PDB_read import PDB_read as pdbr
import copy

class BuildDictionary:
	ala = pdbr('amino_acids/alanine.pdb')
	arg = pdbr('amino_acids/arginine.pdb')
	asn = pdbr('amino_acids/asparagine.pdb')
	asp = pdbr('amino_acids/aspartic_acid.pdb')
	cys = pdbr('amino_acids/cysteine.pdb')
	glu = pdbr('amino_acids/glutamic_acid.pdb')
	gln = pdbr('amino_acids/glutamine.pdb')
	gly = pdbr('amino_acids/glycine.pdb')
	his = pdbr('amino_acids/histidine.pdb')
	ile = pdbr('amino_acids/isoleucine.pdb')
	leu = pdbr('amino_acids/leucine.pdb')
	lys = pdbr('amino_acids/lysine.pdb')
	met = pdbr('amino_acids/methionine.pdb')
	phe = pdbr('amino_acids/phenalalanine.pdb')
	pro = pdbr('amino_acids/proline.pdb')
	ser = pdbr('amino_acids/serine.pdb')
	thr = pdbr('amino_acids/threonine.pdb')
	trp = pdbr('amino_acids/tryptophan.pdb')
	tyr = pdbr('amino_acids/tyrosine.pdb')
	val = pdbr('amino_acids/valine.pdb')
	dic = None

	def __init__(self):
		self.dic = {"A":self.ala.getAtom(),"R":self.arg.getAtom(),"N":self.asn.getAtom(),"D":self.asp.getAtom(),"C":self.cys.getAtom(),"E":self.glu.getAtom(),
					"Q":self.gln.getAtom(),"G":self.gly.getAtom(),"H":self.his.getAtom(),"I":self.ile.getAtom(),"L":self.leu.getAtom(),"K":self.lys.getAtom(),
					"M":self.met.getAtom(),"F":self.phe.getAtom(),"P":self.pro.getAtom(),"S":self.ser.getAtom(),"T":self.thr.getAtom(),"W":self.trp.getAtom(),
					"Y":self.tyr.getAtom(),"V":self.val.getAtom()} #Creates a dictionary with the letter representation of amino acids as key and its atomic structure as values

	def getDic(self):
		return copy.deepcopy(self.dic)