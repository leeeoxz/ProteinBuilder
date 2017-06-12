from ProteinBuilder import ProteinBuilder as buildProt

seq = "VSCEDCPEHCSTQKAQAKCDNDKCVCEPI"

builder = buildProt(seq)

structure = builder.getStructure()

file = open(seq+".pdb","w") #raw_input("Output file (without .pdb): ")

for aa in structure:
	for atom in aa:
		line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(atom[0], atom[1], atom[2], " ", atom[3], " ", atom[4], " ", atom[5], atom[6], atom[7], atom[8], atom[9], " ", " ")
		file.write(line)
		file.write("\n")
file.write("TER")
file.close()