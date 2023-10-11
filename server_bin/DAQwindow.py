import os
import numpy as np
import csv
import string
import argparse

def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("INP",type=str, help="(in PDB format)")
    parser.add_argument("Lwindow", type=int, help="(in PDB format)")
    parser.add_argument('--OutPath', type=str, dest='OutPath', default='./', help='Final Outpt PDB file')
    args = parser.parse_args()
    return args


def get_resscore(filename,window,chain_id):
    p={}
    with open(filename) as result:
        lines = [l for l in result]
        for i in range(len(lines)):
            l = lines[i]
            if (l.startswith('ATOM')) and l[21]==chain_id:
                #only CA
                if l[13:15] != 'CA':
                    continue
                resn=l[22:26].replace(' ','')
                #print('|'+l[30:38]+'|'+l[38:46]+'|'+l[46:54]+'|',resn)
                res_name = l[17:20]
                x=float(l[30:38])
                y=float(l[38:46])
                z=float(l[46:54])
                cd=([x,y,z])
                sco=float(l[60:66])
                l2 = lines[i+1]
                ca_sco = float(l2.split()[2])
                
                p[resn]=[cd,sco,res_name,ca_sco]
                #print('Res',resn,p[resn])
        #Window Score
        for resn in p:
            cnt=0
            sco=0 #p[resn][1] DAQ(AA)
            ca_sco=0 #p[resn][3] DAQ(Ca)
            
            for check_c in range(int(resn)-window,int(resn)+window+1):
                c=str(check_c)
                if c in p:
                    sco=sco+p[c][1]
                    ca_sco=ca_sco+p[c][3]
                    cnt=cnt+1
            p[resn].append(sco/float(cnt)) #p[resn][4]
            p[resn].append(ca_sco/float(cnt)) #p[resn][5]
    return p

def get_pdb(filename):
    p={}
    with open(filename) as result:
        for l in result:
            if(l.startswith('ATOM') and l[13:16] == 'CA '):
                resn=l[22:26].replace(' ','')
                #print('|'+l[30:38]+'|'+l[38:46]+'|'+l[46:54]+'|',resn)
                x=float(l[30:38])
                y=float(l[38:46])
                z=float(l[46:55])
                cd=[x,y,z]
                sco=float(l[60:66])
                AA = l[17:20]
                p[resn]=[cd,AA]
    return p
from collections import defaultdict
def read_pdb_info(filename,chain_id):
    #read each residues for all other informations
    residue_dict=defaultdict(list)
    with open(filename) as result:
        for l in result:
            if l.startswith('ATOM') and l[21]==chain_id:
                chain_name = l[21]
                atom_name = l[12:16]
                x=float(l[30:38])
                y=float(l[38:46])
                z=float(l[46:55])
                resn=l[22:26].replace(' ','')
                residue_dict[resn].append([chain_name,atom_name,x,y,z])
    return residue_dict
def read_chain_set(filename):
    #read each residues for all other informations
    chain_set = set()
    with open(filename) as result:
        for l in result:
            if l.startswith('ATOM'):
                chain_name = l[21]
                chain_set.add(chain_name)
    return chain_set


def save_pdb_with_score(p,residue_dict,filename):
    
    
    output = open(filename, 'w')
    Natm=1
    for resn in p:
            
        sco = p[resn][4]
        ca_sco = p[resn][5]
        current_residue = residue_dict[resn]
        print(current_residue)
        for item in current_residue:
            line='ATOM{:7d} {:4} {:3} {:1}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}\n'.\
                format(Natm,item[1],p[resn][2] ,item[0],int(resn),item[2],item[3],item[4],sco)
            Natm = Natm+1
            output.write(line)
    output.close()

def save_pdb_with_score_filter(p,residue_dict,filename,aa_cut,ca_cut):
    
    
    output = open(filename, 'w')
    Natm=1
    for resn in p:
            
        sco = p[resn][4]
        ca_sco = p[resn][5]
        current_residue = residue_dict[resn]
        print(current_residue,f'AA_SCO= {sco} CA_SCO= {ca_sco}')

        #remove ca_sco < -0.5
        if ca_sco < ca_cut:
            continue
        aa_name = p[resn][2]
        if sco < aa_cut:
            aa_name = 'UNK'
        for item in current_residue:
            line='ATOM{:7d} {:4} {:3} {:1}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}\n'.\
                format(Natm,item[1],aa_name ,item[0],int(resn),item[2],item[3],item[4],sco)
            Natm = Natm+1
            output.write(line)
    output.close()

#=========================================================
def main():
    args = get_args()
    
    # Instantiate the data problem.
    filename1 = args.INP
    raw_score_save_path = args.INP
    pdb_path = raw_score_save_path
    save_path=args.OutPath
    window_size = args.Lwindow

    score_save_path = os.path.join(save_path+"daq_score_w"+str(window_size)+".pdb")
    chain_list = read_chain_set(pdb_path)

    print("%s total different chains :",chain_list)
    
    for chain_name in chain_list:
        score_chain_save_path = os.path.join(save_path+"daq_score_w"+str(window_size)+"_"+str(chain_name)+".pdb")
        score_dict = get_resscore(raw_score_save_path,window_size,chain_name)
        residue_dict = read_pdb_info(pdb_path,chain_name)
        save_pdb_with_score(score_dict, residue_dict,score_chain_save_path)
        #concat all chain visualization together
        #UNK and remove low DAQ(CA) positions
        score_chain_save_path = os.path.join(save_path+"daq_score_w"+str(window_size)+"_"+str(chain_name)+"_filter.pdb")
        score_dict = get_resscore(raw_score_save_path,window_size,chain_name)
        residue_dict = read_pdb_info(pdb_path,chain_name)
        save_pdb_with_score_filter(score_dict, residue_dict,score_chain_save_path,-0.5,-0.5) #AA, CA

    with open(score_save_path,'w') as wfile:
        #Original DAQ(AA)
        wfile.write('MODEL  1\n')
        for chain_name in chain_list:
            score_chain_save_path = os.path.join(save_path+"daq_score_w"+str(window_size)+"_"+str(chain_name)+".pdb")
            with open(score_chain_save_path,'r') as rfile:
                line = rfile.readline()
                while line:
                    wfile.write(line)
                    line = rfile.readline()
        wfile.write('ENDMDL\n')
        #Filtered Model
        wfile.write('MODEL  2\n')
        for chain_name in chain_list:
            score_chain_save_path = os.path.join(save_path+"daq_score_w"+str(window_size)+"_"+str(chain_name)+"_filter.pdb")
            with open(score_chain_save_path,'r') as rfile:
                line = rfile.readline()
                while line:
                    wfile.write(line)
                    line = rfile.readline()
        wfile.write('ENDMDL\n')
            






if __name__ == '__main__':
    main()
