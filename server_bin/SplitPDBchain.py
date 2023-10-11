import os

import sys
import string
import argparse
import subprocess

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

parser = PDBParser()
io = PDBIO()



def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("PDB", type=str, help="INPUT PDB formatfile")
    parser.add_argument('--OutPath', type=str, dest='OutPath', default='./', help='Final Outpt Path for MRC and PDB files')
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    InputPDBpath = args.PDB
    #reso = args.Reso



    if True:
        outpath = args.OutPath
        structure = parser.get_structure("1inp",InputPDBpath )
        pdb_chains = structure.get_chains()
        for chain in pdb_chains:
            chid = chain.get_id()
            #outmrcfile = outpath + 'R'+str(reso)+'_'+chid+'.mrc'
            outpdbfile = outpath + '_'+chid+'.pdb'

            if not os.path.isfile(outpdbfile):
                io.set_structure(chain)
                io.save(outpdbfile)
                #cmd = ['/apps/eman2/EMAN2_2.12/bin/e2pdb2mrc.py','-R',reso,outpdbfile,outmrcfile]
                #res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')


    #except:
    #    print('Failed')



if __name__ == '__main__':
    main()

