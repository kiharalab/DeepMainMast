import os
import subprocess
import numpy as np
from Bio.PDB import *
import random
import string
import argparse
import copy


def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("OUT", type=str, help="Vesper output file")
    parser.add_argument("MODEL", type=str, help="Model(PDB format)")
    parser.add_argument('--OutPath', type=str, dest='OutPath', default='./', help='Out put path. Please use Dir name such as ./DIR/Prob05_ then you will get ./DIR/Prob05_FIT_MODEL*.pdb')
    args = parser.parse_args()
    return args

def main():
    #cmd=['ln','-s',map_path,input_map]
    #res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')

    args = get_args()
    MODEL=args.MODEL
    OutPath=args.OutPath
    if not os.path.isfile(args.OUT):
        print('Can not find',args.OUT)
        return 0
    TopMtx,Trans=load_vesper_out(args.OUT)
    
    model = PDBParser().get_structure("model",MODEL)[0]
    
    print(Trans)
    Num=0
    for rota,tra in zip(TopMtx,Trans):
        rmodel=copy.deepcopy(model)
        print(rota,tra)
        tmtx = np.array(tra).T
        #tmtx = np.array([0.0,0.0,0.0])
        rmtx = np.array([rota[0:3],rota[3:6],rota[6:9]]).T
        
        print(rmtx)
        
        for atm in rmodel.get_atoms():
            atm.transform(rmtx,tmtx)
            
        
        
        outfile=OutPath+'FIT_MODEL'+str(Num)+'.pdb'
        io = PDBIO()
        io.set_structure(rmodel)
        io.save(outfile)
        Num=Num+1
    
def load_vesper_out(file):
    MTX=[]
    TR=[]
    with open(file) as f:
        for l in f:
            #print(l)
            if 'Q=' in l and 'zsco=' in l:
                tmp=[float(i.split('=')[-1].split('{')[-1]) for i in l.replace('}','').replace('{','').split()[5:14]]
                #print(tmp)
                MTX.append(tmp)
                
                tmp=[float(i.split('=')[-1].split('{')[-1]) for i in l.replace('}','').replace('{','').split()[14:17]]
                TR.append(tmp)
    return MTX,TR



    
if __name__ == '__main__':
    main()