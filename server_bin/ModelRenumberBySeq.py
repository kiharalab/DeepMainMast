import os
import numpy as np
import csv
import string
import argparse
from Bio import Align
from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO
from Bio.SeqUtils import seq1
import numpy as np
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB import Model



def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("INP",type=str, help="(in PDB format)")
    parser.add_argument("SEQ", type=str, help="(in FASTA format)")
    parser.add_argument('--OutPath', type=str, dest='OutPath', default='./renumber.pdb', help='Final Outpt PDB file')
    parser.add_argument('--TempFile', type=str, dest='Tmp', default='./tmp.pdb', help='Tmp PDB file')
    args = parser.parse_args()
    return args

def comp_seq_id(ali):
    
    iden_pos = sum(res1 == res2 for res1,res2 in zip(ali[0],ali[2]))
    total_pos = sum(1 for res1,res2 in zip(ali[0],ali[2]) if res1 !='-')
    if total_pos == 0:
        return 0
    else:
        return iden_pos/total_pos


# Dictionary to map three-letter code to one-letter code
three_to_one = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G', 
                'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 
                'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}

chain_table = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

def edit_pdb_file(input_file, output_file):
    chain_counter = 1
    with open(input_file, 'r') as inp, open(output_file, 'w') as out:
        for line in inp:
            if line.startswith('ATOM'):
                chain_id = line[21]
                #new_chain_id = f'{chain_id}{chain_counter}'
                new_chain_id = chain_table[chain_counter-1]
                new_line = line[:21] + new_chain_id + line[22:]
                out.write(new_line)
            elif line.startswith('TER'):
                chain_counter += 1
                out.write(line)
            elif line.startswith('END'):
                continue
            else:
                out.write(line)


#=========================================================
def main():
    args = get_args()
    
    # Instantiate the data problem.
    init_path = args.INP
    tmp_path = args.Tmp
    #edit by TER
    edit_pdb_file(init_path,tmp_path)


    seq_path = args.SEQ
    save_path=args.OutPath
    sequences = []
    cid = 0
    seq_tbl ={}
    chain_list = []
    start_tbl={}
    start_tbl['A']=1
    with open(seq_path) as handle:
        for record in SeqIO.parse(handle,"fasta"):
            sequence = record.seq
            chain_name = chain_table[cid]
            seq_tbl[chain_name] = sequence
            cid = cid + 1
            next_chain_name = chain_table[cid]
            start_tbl[next_chain_name]=start_tbl[chain_name]+len(sequence)

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("structure",tmp_path)
    print(len(structure[0]))
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score=1
    aligner.mismatch_score=0
    aligner.target_internal_open_gap_score=-1


    final_str = Model.Model(1)

    #multiple chains in the sequence and model
    for chain_name in seq_tbl:
        print(chain_name,seq_tbl[chain_name])
        seq = seq_tbl[chain_name]
        for chain in structure[0]:
            pdb_sequence = ''
            residue_mapping = {}
            for residue in chain:
                res_name=residue.get_resname()
                if res_name in three_to_one:  # Check if residue name is in dictionary
                    pdb_sequence += three_to_one[res_name]
                    residue_mapping[len(pdb_sequence)] = residue

            alignments = aligner.align(seq,pdb_sequence)
            #print(pdb_sequence)
            #print(alignments[0])
            ali = str(alignments[0]).split('\n')
            print(ali)
            #seqID=comp_seq_id(ali)
            seqID=alignments[0].score/len(pdb_sequence)
            print(f'SeqID= {seqID}')
            if seqID != 1.0:
                continue
            #create new residue numbers
            new_resnum = {}
            new_number = start_tbl[chain_name]
            for i in range(len(ali[0])):
                if ali[0][i] == ali[2][i]:
                    print(ali[0][i],ali[2][i],new_number)
                    new_resnum[residue_mapping[i+1].id[1]] = new_number
                    new_number += 1
            for residue in new_resnum:
                #print(residue_mapping[residue])
                residue_mapping[residue].id = (' ',new_resnum[residue],' ')
            final_str.add(chain)
    io = PDBIO()
    io.set_structure(final_str)
    io.save(save_path)






if __name__ == '__main__':
    main()
