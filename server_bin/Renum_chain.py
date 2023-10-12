import os
import copy
import sys
import string
import argparse
import subprocess

from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio import PDB, SeqIO,AlignIO
from Bio.Align import PairwiseAligner

def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("FASTA", type=str, help="INPUT PDB formatfile")
    parser.add_argument("PDB", type=str, help="INPUT PDB formatfile")
    parser.add_argument('--OutPath', type=str, dest='OutPath', default='./renumbered.pdb', help='Outpt Path of the renumbered PDB file')
    args = parser.parse_args()
    return args

def split_structures(filename,outpath):

    #Split file
    Nmodel=0
    sublines=''
    with open(filename) as f:
        lines=[ i for i in f]
        for l in lines:
            if l.startswith('TER'):
                #Save tmp model
                if len(sublines)>0:
                    with open(f'{outpath}_{Nmodel}','w') as OUT:
                        OUT.write(sublines)
                        print(f'writing...{outpath}_{Nmodel}')
                        Nmodel = Nmodel + 1
                sublines=''
            else:
                if l.startswith('ATOM'):
                    sublines = sublines + l
        if len(sublines)>0:
            with open(f'{outpath}_{Nmodel}','w') as OUT:
                OUT.write(sublines)
                print(f'writing...{outpath}_{Nmodel}')
                Nmodel = Nmodel + 1

    parser = PDB.PDBParser(QUIET=True)
    structures=[]
    for i in range(Nmodel):
        structure = parser.get_structure(f'model{i}',f'{outpath}_{i}')
        structures.append(structure)
    return structures


def load_multiple_seq(filename):
    results=[]
    seq_data=''
    with open(filename) as f:
        lines=[ i for i in f]
        for l in lines:
            if l.startswith('>'):
                if len(seq_data)>0:
                    results.append(seq_data)
                seq_data=''
            else:
                seq_data = seq_data + l.strip().replace(' ','')
        if len(seq_data)>0:
            results.append(seq_data)
    return results

def get_sequence_from_model(structure):
    seq_data=''
    resnum_tbl=[]
    for chain in structure:
        for residue in chain:
            if PDB.is_aa(residue,standard=True):
                seq_data += PDB.Polypeptide.three_to_one(residue.get_resname())
                resnum_tbl.append(residue.id[1])
    return seq_data,resnum_tbl

def renumber_pdb(ali,model,start_num,res_tbl,cid):
    new_structure = PDB.Structure.Structure('new_model')
    new_model = PDB.Model.Model(0)
    new_chain = PDB.Chain.Chain(cid)
    lines = str(ali).split('\n') #Bug in Biopython??
    ali1 = lines[0]
    ali2 = lines[2]
    fasta_index = start_num
    pdb_index = 0
    fasta_tbl=[]
    pdb_tbl=[]
    for i in range(len(ali1)):
        if ali1[i] == ali2[i]:
            fasta_tbl.append(fasta_index)
            pdb_tbl.append(res_tbl[pdb_index])
            #pdb_tbl.append(pdb_index)
        if ali1[i] != '-':
            fasta_index += 1
        if ali2[i] != '-':
            pdb_index += 1
    #print(pdb_tbl,fasta_tbl)
    #print(model)
    for res_num,new_res_num in zip(pdb_tbl,fasta_tbl):
        #print('**',res_num,new_res_num)
        for chain in model[0]:
            if res_num in chain:
                residue = chain[res_num]
                #create new residue
                new_residue = PDB.Residue.Residue(
                    (' ',new_res_num, ' '),
                    residue.resname,
                    residue.segid
                )
                #Copy Atoms
                for atom in residue:
                    new_residue.add(atom.copy())
                new_chain.add(new_residue)
    new_model.add(new_chain)
    new_structure.add(new_model)
    return new_structure
        

def main():
    args = get_args()
    InputPDBpath = args.PDB
    InputFASTApath = args.FASTA
    OutPath = args.OutPath
    ChainIdx=string.ascii_uppercase + string.ascii_lowercase + string.ascii_uppercase + string.ascii_lowercase #Max 24*2

    io = PDBIO()
    if not os.path.isfile(InputFASTApath):
        print(f'Missing FASTA file: {InputFASTApath}')
        return True
    
    if not os.path.isfile(InputPDBpath):
        print(f'Missing PDB file: {InputPDBpath}')
        return True

    if os.path.isfile(OutPath):
        print(f'Already Existing the file: {OutPath}')
        return True


    structures = split_structures(InputPDBpath,OutPath)

    pdb_sequences = []
    rtbls=[]
    for model in structures:
        seq_data,rtbl = get_sequence_from_model(model[0])
        pdb_sequences.append(seq_data)
        rtbls.append(rtbl)

    fasta_sequences = load_multiple_seq(InputFASTApath)
    start_tbl=[1]
    for i in range(len(fasta_sequences)):
        start_tbl.append(len(fasta_sequences[i])+start_tbl[i])
    print('START:',start_tbl)

    print('FASTA',fasta_sequences)
    print('PDB:',pdb_sequences)

    #Make alignments
    aligner = PairwiseAligner()


    merged_structure = PDB.Structure.Structure('merged')
    model_number=0
    for i in range(len(fasta_sequences)):
        ch=ChainIdx[i]
        for j in range(len(pdb_sequences)):
            fseq = fasta_sequences[i]
            pseq = pdb_sequences[j]

            alignments = aligner.align(fseq,pseq)
            print(f'FASTA{i}:PDB{j}',alignments[0].score,len(pseq))
            if alignments[0].score == len(pseq):
                print(f'Assign chain {ch}')
                #renumber residue IDs
                alignment = alignments[0]
                new_data = renumber_pdb(alignment,structures[j],start_tbl[i],rtbls[j],ch)

                new_model = PDB.Model.Model(model_number)
                for tmp_model in new_data:
                    for chain in tmp_model:
                        new_model.add(chain.copy())
                merged_structure.add(new_model)
                model_number += 1
                break
    io = PDB.PDBIO()
    io.set_structure(merged_structure)
    io.save(OutPath)

if __name__ == '__main__':
    main()