from PDB_read import PDB_read
from SA import SimAnne
from RMSD import rmsd
import math
import numpy as np
import copy
import random

r= rmsd("1plx.pdb","YGGFM.pdb","backbone")
r.rotine()