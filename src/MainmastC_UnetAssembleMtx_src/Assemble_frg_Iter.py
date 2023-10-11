#!/usr/bin/env python
# coding: utf-8

# In[158]:


import os
import subprocess
import numpy as np
#print('BIO')
import Bio
#import networkx as nx
import random
#import matplotlib.pyplot as plt
#print('PULP')
from pulp import *
from pulp.apis import PULP_CBC_CMD
from ortools.linear_solver import pywraplp
from ortools.init import pywrapinit
import string
import argparse
from ortools.sat.python import cp_model
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
from numba import jit


# In[199]:

def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("LIST", type=str, help="Output from MainmastC_fastfragMtx")
    parser.add_argument('--Nstock', type=int, dest='Nstock', default=3000, help='Number of Stocks for SAT solver')
    parser.add_argument('--SecAssemble', type=int, dest='SecCp', default=300, help='Computatinal Time Limit for Assemble')
    parser.add_argument('--OutPath', type=str, dest='OutNamePath', default='final.pdb', help='Final Outpt PDB file')
    args = parser.parse_args()
    return args

def create_data_model(filename):
    """Stores the data for the problem."""
    data = {}
    data['fragment']=[]
    data['sequence']=[]
    data['connection']=[]
    data['node']=[]
    seq=[]
    with open(filename) as f:
        for l in f:
            if l.startswith('FRG'):
                fid   = int(l.split()[1]) 
                rawsco = int(float(l.split()[3])*100) ##Int
                zsco = int(float(l.split()[5])*100) ##Int
                ali = l.split()[9]
                data['fragment'].append([fid,rawsco,zsco,ali])
                
            if l.startswith('SEQ'):
                x = int(l.split()[1])
                seq.append([x,l.split()[3],int(l.split()[2])]) #number, AA, chain
                
            if l.startswith('MTX'):                
                x = int(l.split()[1])
                if len(l.split()) == 2:
                    y = '-1,-1'
                else:
                    y = l.split()[2]
                data['connection'].append([x,y])
            if l.startswith('NODE'):
                x = float(l.split()[1])
                y = float(l.split()[2])
                z = float(l.split()[3])
                data['node'].append([x,y,z])
        
        Lseq = seq[-1][0] #Length of  sequence

        print('Lseq=',Lseq)
        data['sequence'] = [ [0,'XXX',-1] for i in range(Lseq+1)]
    
        for row in seq:
            data['sequence'][row[0]]=[row[0],row[1],row[2]]
    
    return data

def MakeSubset(fdata,alldata,cdata,Nsub):
    import random
            
    Total = len(fdata)
    Nround = int(Total/Nsub) + 1
    out=[]
    labels = [i[0] for i in fdata] #ID
    OUT = [[] for i in range(Nround)]
    random.shuffle(labels)
    for i in range(Nround):        
        start = i * Nsub
        end = start + Nsub
        print('###Subset:',i, start, end)
        out.append([labels[j] for j in range(start,end) if j < Total])
        tbl = out[-1]
        #print(tbl)
        #update fragment tbl
        for j in tbl:
            #print(i,j,cdata[j])
            
            tabu = [int(k) for k in cdata[j][1].split(',')[:-1]]
            tabu2 = set(tbl) & set(tabu) #common
            #print(len(tabu2),len(tabu),len(tabu))
            OUT[i].append(alldata[j]+list(tabu2))
    return OUT,Nround

def Assemble(d,sec):
    model = cp_model.CpModel()

    # Variables
    #sequence - position
    Nf = len(d)
    x = []
    
    idtbl = [i[0] for i in d] #order -> ID
    id2tbl =[-1 for i in range(max(idtbl)+1)] #ID-> Order
    for i in range(len(idtbl)):
        id2tbl[idtbl[i]]=i

    #[fid,rawsco,zsco,ali,tabu]
    for i in range(len(idtbl)):
        x.append(model.NewBoolVar(f'x[{i}]'))
        #print(i[4:10])
                
    #Const
    print('#Adding constraints')
    for i in range(len(idtbl)):
        #print(i,d[i][4:10])
        model.Add(sum(x[id2tbl[j]] for j in d[i][4:] if j>=0)==0).OnlyEnforceIf(x[i])
        
        
    #return True
    print('#Setting Objective Score')
    # Objective
    objective_terms = []
    #[fid,rawsco,zsco,ali,tabu]
    for i in range(Nf):
        objective_terms.append(d[i][1] * x[i]) #Sum of Raw Scores
    model.Maximize(sum(objective_terms))

    print('#Start Solving...')
    # Solve
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = sec
    solver.parameters.num_search_workers = 8
    solver.parameters.log_search_progress = True
    solution_printer = cp_model.ObjectiveSolutionPrinter()
    #status = solver.Solve(model)
    status = solver.SolveWithSolutionCallback(model, solution_printer)
    results=[]
    if status == pywraplp.Solver.OPTIMAL or pywraplp.Solver.FEASIBLE:
        for i in range(len(idtbl)):
            if solver.BooleanValue(x[i]):
                #print(i,'ID=',idtbl[i]) #order->ID
                results.append(d[i])
    return results

def ShowPDB(results,d,filename):
    chainIdx=string.ascii_uppercase + string.ascii_lowercase + string.ascii_uppercase + string.ascii_lowercase #Max 24*2
    outlines=''
    Natm = 1
    #Sequence Order
    for pos in d['sequence']: #number, AA, chain
        print(pos)
        
        #Check Coordinates
        x=np.array([])
        y=np.array([])
        z=np.array([])
        for frg in results:
            #print(frg)
            ali = [int(i) for i in frg[3].split(',')[0:-1]]
            #print(pos[0], ali)
            for ca in range(int(len(ali)/2)):
                if ali[2*ca+1] == pos[0]:
                    #print(ali[ca*2],ali[2*ca+1]) #node, sequence
                    x=np.append(x,d['node'][ali[ca*2]][0])
                    y=np.append(y,d['node'][ali[ca*2]][1])
                    z=np.append(z,d['node'][ali[ca*2]][2])
        if len(x)==0: #Skip
            continue
        x_ave = np.average(x)
        y_ave = np.average(y)
        z_ave = np.average(z)
        #print(x, x_ave)
        
        chainID = chainIdx[pos[2]]
        outlines +='ATOM{:7d}  CA  {:3} {:1}{:4d}'.format(Natm,pos[1],chainID,pos[0])
        outlines +='    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}\n'.format(x_ave,y_ave,z_ave,1.0)
        Natm = Natm+1
    
    with open(filename,'w') as out:
        out.write(outlines)

        
        

#=========================================================
def main():
    args = get_args()
    """Entry point of the program."""
    # Instantiate the data problem.
    filename = args.LIST
    #filename = 'tmp2.txt'
    data = create_data_model(filename)
    print('#Nfragments=',len(data['fragment']))
    print('#Nnode=',len(data['node']))

    #Assemble Part
    Nsub = args.Nstock
    Sec = args.SecCp
    TMP_FRAG = data['fragment']
    for it in range(20):
        #Make Subset
        print('Itre',it)
        SubSet, Nrnd = MakeSubset(TMP_FRAG,data['fragment'],data['connection'],Nsub)
        print(it,'Nround=',Nrnd)
        TMP_FRAG=[]
        for sub in range(Nrnd):
            print('##Sub',sub+1,'/',Nrnd)
            asb_results=Assemble(SubSet[sub],Sec)
            
            #show
            for i in asb_results:
                TMP_FRAG.append(i)
        print('#DONE',it,'Nfrag=',len(TMP_FRAG))
        if Nrnd == 1:
            print('Finished')
            break
    #Make Final Model
    outputfile = args.OutNamePath
    ShowPDB(TMP_FRAG,data,outputfile)
    
    
if __name__ == '__main__':
    main()
