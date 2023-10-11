from TEMPy.protein.structure_blurrer import StructureBlurrer
from TEMPy.protein.structure_parser import PDBParser
import sys

if len(sys.argv) != 4:
	print("Incorrect arguments")
	exit(0)

reso = float(sys.argv[1])
pdb_file = sys.argv[2]
mrc_file = sys.argv[3]

structure_instance = PDBParser.read_PDB_file("X", pdb_file, hetatm=False, water=False)
blurrer = StructureBlurrer()
modelmap = blurrer.gaussian_blur_real_space(structure_instance, reso, sigma_coeff=0.187, densMap=False, normalise=True) 
modelmap.write_to_MRC_file(mrc_file)