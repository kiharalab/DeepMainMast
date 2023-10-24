import os
import subprocess
import numpy as np
import argparse
#Format check

def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("SEQ", type=str, help="Input Sequence file")
    parser.add_argument('--OutPath', type=str, dest='OutPath', default='out.fasta', help='Outpt FASTA file')
    args = parser.parse_args()
    return args

AA20='ARNDCQEGHILKMFPSTWYV'

def main():
    args = get_args()
    inputpath = args.SEQ
    outputpath = args.OutPath
    outlines=''
    seq=''
    with open(inputpath) as f:
        lines = [l for l in f]
        #print(lines)
        for l in lines:
            l = l.strip()
            #Header lines
            if l.startswith('>'):
                if len(seq)>0:
                    outlines += f'{seq}\n'
                outlines += f'{l}\n'
                seq=''
            else:
                #check 20AA or not
                for aa in l:
                    if aa in AA20:
                        seq+=aa
                        #print(aa)
                    else:
                        print(f'Unknown AA type {aa}')
        if len(seq)>0:
            outlines += f'{seq}\n'
        print(outlines)
    with open(outputpath,'w') as out:
        out.write(outlines)






if __name__ == '__main__':
    main()
