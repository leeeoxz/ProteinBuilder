from PDB_read import PDB_read as pdbr
from Rotational import getAngle as angle
import numpy as np
import matplotlib.pyplot as plt
import os

f = pdbr("1eny.pdb")
ca = f.getAlpha() #gets all Alpha Carbon
n = f.getN() #Gets the N from backbone
c = f.getC() #Gets the C from backbone
aa = f.getAA() # gets all AA names
phis = []
psis = []

f = open("1eny.txt","w")
f.write("Amino"+"\t"+"Phi"+"\t"+"Psi\n")

def isNumber(atom):
		try:
			int(atom[4])
			return True

		except ValueError:
			return False

if isNumber(ca[0]):
	start = 5
else:
	start = 6


for x,item in enumerate(aa):
	if x == 0:
		phi = 360.00
		psi = angle(n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3],n[x+1][start:start+3])
		apsi = psi.getAngle()
		psis.append(apsi)
		phis.append(phi)
		f.write("{:5s} {:7.2f} {:7.2f}".format(item,phi,apsi)+"\n")
	elif x == len(aa) - 1:
		phi = angle(c[x-1][start:start+3],n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3])
		psi = 360.00
		aphi = phi.getAngle()
		psis.append(psi)
		phis.append(aphi)
		f.write("{:5s} {:7.2f} {:7.2f}".format(item,aphi,psi))
	else:
		phi = angle(c[x-1][start:start+3],n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3])
		psi = angle(n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3],n[x+1][start:start+3])
		apsi = psi.getAngle()
		aphi = phi.getAngle()
		psis.append(apsi)
		phis.append(aphi)
		f.write("{:5s} {:7.2f} {:7.2f}".format(item,aphi,apsi)+"\n")
plt.plot(phis, psis, 'ro', color = 'grey', ms = 5.0)
plt.arrow( -180, 0, 360, 0 )   # Creates an arrow
plt.arrow( 0, -180, 0, 360 )   # Creates an arrow
plt.ylabel('PSI')
plt.xlabel('PHI')
plt.title('1eny.pdb Ramachandran')
plt.axis([-180, 180, -180, 180])
plt.savefig("1eny.eps", dpi=1.000, facecolor='w', edgecolor='w',
        orientation='portrait', papertype="letter", format="eps",
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)