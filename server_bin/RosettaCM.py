import os
import subprocess
import numpy as np
from Bio.PDB import *
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import PDBParser

import random
import string
import argparse
import copy

parser = PDBParser()
io = PDBIO()

def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("FASTA", type=str, help="Target Sequence")
    parser.add_argument("MODEL", type=str, help="Model(PDB format)")
    parser.add_argument("MAP", type=str, help="MAP")
    parser.add_argument('--OutPath', type=str, dest='OutPath', default='./OutPutDir', help='Out put path. Please use Dir name such as ./OUTPUTs')
    parser.add_argument('--XMLPath', type=str, dest='XMLPath', default='/home/users/', help='Copy *.xml and *.sh files')
    parser.add_argument('--PulchraPath', type=str, dest='PulchraPath', default='/home/users/bin/pulchra', help='Path of PULCHRA exe file')
    args = parser.parse_args()
    return args

def main():
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
             'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    #cmd=['ln','-s',map_path,input_map]
    #res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')

    args = get_args()
    MODEL=args.MODEL
    OutPath=args.OutPath
    FASTA=args.FASTA
    MAP=args.MAP
    PulchraPath=args.PulchraPath
    if not os.path.isfile(FASTA):
        print('Can not find',FASTA)
        return 0
    if not os.path.isfile(MAP):
        print('Can not find',MAP)
        return 0
    if not os.path.isfile(MAP):
        print('Can not find',MAP)
        return 0
    
    if os.mkdir(OutPath):
        print('Directory is existing..',OutPath)
        return 0
    
    
    model = PDBParser().get_structure("model",MODEL)[0]
    
    id_name,seq = ReadFasta(FASTA)
    
    print(id_name,seq)
    #Multi-chain
    #get first chain
    #for chain in model:
    #    break
    temp_seq=''
    SeqPos=0
    Residues = model.get_residues()
    print(Residues)
    for pos in range(len(seq)):
        if seq[pos] != '/':#Chain Break
            SeqPos = SeqPos + 1
        else:
            temp_seq=temp_seq+'/'
            continue

        print(SeqPos,seq[pos])
        found=False
        for chain in model:
            if SeqPos in chain:
                print('Detected:',d3to1[chain[SeqPos].resname])
                temp_seq=temp_seq+d3to1[chain[SeqPos].resname]
                found=True
        if found == False:
            temp_seq=temp_seq+'-'
            
    #Make Alignment file
    alignment_file=OutPath+'/alignment.txt'
    lines='## 1XXX_ 1tmpA_thread\n'
    lines=lines+'# hhsearch\n'
    lines=lines+'scores_from_program: 0 1.00\n'
    #remove chain-break in the alignment
    seq_ali = seq.replace('/','')
    temp_seq_ali = temp_seq.replace('/','')
    lines=lines+'0 '+seq_ali+'\n'
    lines=lines+'0 '+temp_seq_ali+'\n--\n'
    #lines=lines+'0 '+seq+'\n'
    #lines=lines+'0 '+temp_seq+'\n--\n'
    print(seq)
    print(temp_seq)
    with open(alignment_file,'w') as out:
        out.write(lines)
        
        
    out_map = OutPath + '/inputmap.map'
    cmd=['cp',MAP,out_map]
    res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
    
    out_fasta = OutPath + '/seq.fasta'
    #Rosetta Format fasta file
    lines='>Target\n'+seq+'\n'
    with open(out_fasta,'w') as out:
        out.write(lines)
    #cmd=['cp',FASTA,out_fasta]
    #res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
    
    #Prepare template PDB file for each chain
    structure = parser.get_structure("1inp",MODEL)
    pdb_chains = structure.get_chains()
    chain_list =[]
    file_list=[]
    for chain in pdb_chains:
        chid = chain.get_id()
        chain_list.append(chid)
        outpdbfile = OutPath + '/input_'+chid+'.pdb'
        io.set_structure(chain)
        io.save(outpdbfile)
        #cmd=['bin/pulchra','-s',outpdbfile]
        cmd=[PulchraPath,'-s',outpdbfile]
        res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
        fullatm_file = OutPath + '/input_'+chid+'.rebuilt.pdb'
        file_list.append(fullatm_file)

    #Merge *.rebuilt.pdb
    out_model = OutPath + '/1tmpA.pdb'
    with open(out_model,'w') as outf:
        for fname in file_list:
            with open(fname) as infile:
                for line in infile:
                    outf.write(line)

    #out_model = OutPath + '/input.pdb'
    #cmd=['cp',MODEL,out_model]
    #res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
    
    #cmd=['pulchra','-s',out_model]
    #res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
    
    #out_model = OutPath + '/input.rebuilt.pdb'
    #out_model2 = OutPath + '/1tmpA.pdb'
    #cmd=['cp',out_model,out_model2]
    #res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
    
    XMLPath=args.XMLPath
    out_model = XMLPath + '/C_rosettaCM.sh'
    out_model2 = OutPath + '/.'
    cmd=['cp',out_model,out_model2]
    res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
    
    out_model = XMLPath + '/C_rosettaCM.xml'
    out_model2 = OutPath + '/.'
    cmd=['cp',out_model,out_model2]
    res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
    
    out_model = XMLPath + '/A_setup.sh'
    out_model2 = OutPath + '/.'
    cmd=['cp',out_model,out_model2]
    res=subprocess.run(cmd,stdout=subprocess.PIPE,encoding='utf-8')
            

    
def ReadFasta(file):
    seq=''
    name=''
    Nchain = 0
    with open(file) as f:
        for l in f:
            l=l.strip()
            if l.startswith('>'):
                name=l
                Nchain = Nchain + 1
                if Nchain > 1:
                    seq = seq + '/' #chain break
            else:
                seq = seq + l
    return name,seq
           


    
if __name__ == '__main__':
    main()
