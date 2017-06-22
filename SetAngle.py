# Leonardo Santos #
# June 2017 #

from PDB_read import PDB_read as pdbr
from Rotational import getAngle as angle
import numpy as np
# Creates a .txt file with the dihedral angles and plost a Ramachandran's plot #

import matplotlib.pyplot as plt
import os

name = "YGGFM.pdb" #raw_input("File name: ")
f = pdbr(name)
ca = f.getAlpha() #gets all Alpha Carbon
n = f.getN() #Gets the N from backbone
c = f.getC() #Gets the C from backbone
aa = f.getAA() # gets all AA names
phis = [] #phi rotational angles
psis = [] #psi rotational angles

f = open(name.replace("pdb","txt"),"w") #opens a file to return phi and psi angles per amino acid
f.write("Amino"+"\t"+"Phi"+"\t"+"Psi\n") #header Amino   Phi   Psi

def isNumber(atom):
		try: #if havent chain letter
			int(atom[4])
			return True

		except ValueError: #if pdb have the chain letter
			return False

if isNumber(ca[0]): #check if pdb file has the chain letter
	start = 5 #if it doesnt have, the coordinates starts at 5 index
else:
	start = 6 #if it has, the coordinates starts at 6 index


for x,item in enumerate(aa): #for each amino acid and its position
	if x == 0: #if its the first amino acid Phi value is 360.00 by default and calculates only Psi
		phi = 360.00 #default value
		psi = angle(n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3],n[x+1][start:start+3]) #calculate angle
		apsi = psi.getAngle() #return angle
		psis.append(apsi) #includes psi value on a list
		phis.append(phi) #includes phi value on a list
		f.write("{:5s} {:7.2f} {:7.2f}".format(item,phi,apsi)+"\n")
	elif x == len(aa) - 1: #if its the last amino acid Psi value is 360.00 by default and calculates only Phi
		phi = angle(c[x-1][start:start+3],n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3]) #calculate angle
		psi = 360.00 #default value
		aphi = phi.getAngle() #return angle
		psis.append(psi) #includes psi value on a list
		phis.append(aphi) #includes phi value on a list
		f.write("{:5s} {:7.2f} {:7.2f}".format(item,aphi,psi))
	else: #calculate both Psi and Phi angles
		phi = angle(c[x-1][start:start+3],n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3]) #calculate angle
		psi = angle(n[x][start:start+3],ca[x][start:start+3],c[x][start:start+3],n[x+1][start:start+3]) #calculate angle
		apsi = psi.getAngle() #return angle
		aphi = phi.getAngle() #return angle
		psis.append(apsi) #includes psi value on a list
		phis.append(aphi) #includes phi value on a list
		f.write("{:5s} {:7.2f} {:7.2f}".format(item,aphi,apsi)+"\n")
f.close()

#RAMACHANDRAN PLOT#

plt.plot(phis, psis, 'ro', color = 'grey', ms = 5.0)
plt.arrow( -180, 0, 360, 0 )   # Creates an arrow 
plt.arrow( 0, -180, 0, 360 )   # Creates an arrow
plt.ylabel('PSI')
plt.xlabel('PHI')
plt.title(name+' Ramachandran')
plt.axis([-180, 180, -180, 180])
plt.savefig(name.replace("pdb","eps"), dpi=1.000, facecolor='w', edgecolor='w',
        orientation='portrait', papertype="letter", format="eps",
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)