#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import subprocess
import numpy as np
import Bio
import random
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

def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("LIST", type=str, help="Output from MainmastC_fastfragMtx")
    parser.add_argument('--Nstock', type=int, dest='Nstock', default=3000, help='Number of Stocks for SAT solver')
    parser.add_argument('--SecAssemble', type=int, dest='SecCp', default=300, help='Computatinal Time Limit for Assemble')
    parser.add_argument('--OutPath', type=str, dest='OutNamePath', default='final.pdb', help='Final Outpt PDB file')
    parser.add_argument('--Ncpu', type=int, dest='Ncpu', default=8, help='Number of CPU cores')
    parser.add_argument('--Niter', type=int, dest='Niter', default=10, help='Max Number of Iterations')
    parser.add_argument('--HomWeight', type=float, dest='Weight', default=1.0, help='Weight for penalize Inconsistent pairs')
    parser.add_argument('--Keep', type=float, dest='Keep', default=0.9, help='Fraction of Original Fragment')
    args = parser.parse_args()
    return args


def create_data_model(filename):
    """Stores the data for the problem."""
    data = {}
    data['fragment']=[]
    data['sequence']=[]
    data['connect_diff']=[]
    data['connect_same']=[]
    data['pair']=[]
    data['node']=[]
    data['homo']=[]
    seq=[]
    with open(filename) as f:
        for l in f:
            if l.startswith('FRG'):
                mid   = int(l.split()[1]) 
                fid   = int(l.split()[3]) 
                cid   = int(l.split()[5]) 
                rawsco = int(float(l.split()[7])*100) ##Int
                zsco = int(float(l.split()[9])*100) ##Int
                ali = l.split()[13]
                data['fragment'].append([mid,fid,cid,rawsco,zsco,ali])
                
            if l.startswith('SEQ'):
                x = int(l.split()[1])
                seq.append([x,l.split()[3],int(l.split()[2])]) #number, AA, chain
                
            if l.startswith('MTX_SAME'):                
                x = int(l.split()[1])
                if len(l.split()) == 2:
                    y = [-1]
                else:
                    y = [int(i) for i in l.split()[2].split(',') if i != '']
                data['connect_same'].append([x,y])
            if l.startswith('MTX_DIFF'):                
                x = int(l.split()[1])
                if len(l.split()) == 2:
                    y = [-1]
                    #y = '-1,-1'
                else:
                    #y = l.split()[2]
                    y = [int(i) for i in l.split()[2].split(',') if i != '']
                data['connect_diff'].append([x,y])
            if l.startswith('NODE'):
                x = float(l.split()[1])
                y = float(l.split()[2])
                z = float(l.split()[3])
                data['node'].append([x,y,z])
            if l.startswith('PAIR'):#Updated Fid based
                x1 = int(l.split()[1].split(',')[0])
                x2 = int(l.split()[1].split(',')[1])
                y1 = int(l.split()[2].split(',')[0])
                y2 = int(l.split()[2].split(',')[1])
                pen = int(l.split()[3])
                data['pair'].append([x1,x2,y1,y2,pen])
        
        Lseq = seq[-1][0] #Length of  sequence

        print('Lseq=',Lseq)
        data['sequence'] = [ [0,'XXX',-1] for i in range(Lseq+1)]
    
        for row in seq:
            data['sequence'][row[0]]=[row[0],row[1],row[2]]
    
    return data


# In[4]:


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


# In[5]:


def MergeSubset(tmp_results,fdata,alldata,cdata):   
    import random
    merged_list = []
    [merged_list.append(x) for x in tmp_results + fdata if x not in merged_list]
    #merged_list = tmp_results + fdata
    #remove redundant
    
    
    labels = [i[0] for i in merged_list] #ID
    #random.shuffle(labels)
    tbl = labels #selected IDs
    OUT=[]
    for j in tbl:
        tabu = [int(k) for k in cdata[j][1].split(',')[:-1]]
        tabu2 = set(tbl) & set(tabu) #common
        #print(cdata[j][1],tabu[-2],tabu[-1])
        #print(tabu)
        OUT.append(alldata[j]+list(tabu2))
    return OUT


# In[6]:


def Assemble(d,sec,Ncpu,rate):

    Wpen = -1
    model = cp_model.CpModel()

    # Variables
    #[mid,fid,cid,rawsco,zsco,ali]
    tmp_fid=[]
    tmp_cid=[]
    idtbl = {}
    for mid,fid,cid,rawsco,zsco,ali in d['fragment']:
        tmp_fid.append(fid)
        tmp_cid.append(cid)
        idtbl[(fid,cid)]=mid
    max_fid = max(tmp_fid)
    max_cid = max(tmp_cid)

    print('MaxFid= {:d} MaxCid= {:d}'.format(max_fid,max_cid))
    assign = {}

    for f in range(max_fid+1):
        for c in range(max_cid+1):
            assign[(f,c)] = model.NewBoolVar('assign_%i_%i'%(f,c))
            print(assign[(f,c)])
            if not (f,c) in idtbl:
                print("Set False",assign[(f,c)])
                model.Add(assign[(f,c)]==False)
                
                
    #Const1
    print('#Adding constraints')
    for f in range(max_fid+1):
        model.Add(sum(assign[(f,c)] for c in range(max_cid+1))<=1) #Up to one or zero cid
    #Const2 MTX_SAME
    for i,j in d['connect_same']:#fid, others
        for c in range(max_cid+1):
            model.Add(sum(assign[(k,c)] for k in j if k>=0)==0).OnlyEnforceIf(assign[(i,c)])
    
    print('#Adding constraints:Pairs')
    
    print('#Setting Objective Score')
    # Objective
    objective_terms = []
    #[mid,fid,cid,rawsco,zsco,ali]
    SMTX=[[0]*(max_cid+1) for i in range(max_fid+1)]
    for mid,fid,cid,rawsco,zsco,ali in d['fragment']:
        SMTX[fid][cid] = rawsco

    for mid,fid,cid,rawsco,zsco,ali in d['fragment']:
        objective_terms.append(assign[(fid,cid)] * SMTX[fid][cid]) #Sum of Raw Scores

    #PAIRs
    #penalities #x1,x2,y1,y2,pen
    same={}
    for f in range(max_fid+1):
        for c in range(f+1,max_fid+1):
            same[(f,c)] = model.NewBoolVar('same_%i_%i'%(f,c))
    cid={}
    for f in range(max_fid+1):
        cid[(f)] = model.NewIntVar(-1,max_cid+1,'cid_%i'%(f))
    for f in range(max_fid+1):
        for c in range(max_cid+1):
            model.Add(cid[(f)]==c).OnlyEnforceIf(assign[(f,c)])
            model.Add(cid[(f)]!=c).OnlyEnforceIf(assign[(f,c)].Not())
        #model.Add(sum(assign[(f,c)] for c in range(max_cid+1))==0).OnlyEnforceIf(cid[(f)]==-1)
        

    for f1 in range(max_fid+1):
        for f2 in range(f1+1,max_fid+1):
            #same chain
            model.Add(cid[(f1)]==cid[(f2)]).OnlyEnforceIf(same[(f1,f2)])
            model.Add(cid[(f1)]>=0).OnlyEnforceIf(same[(f1,f2)])
            model.Add(cid[(f1)]!=cid[(f2)]).OnlyEnforceIf(same[(f1,f2)].Not())
            #model.Add(sum(assign[(f1,c)] for c in range(max_cid+1))==same[(f1,f2)])
    
    pair={}
    for x1,x2,y1,y2,pen in d['pair']:
        #print(x1,x2,y1,y2)
        pair[(x1,x2,y1,y2)] = model.NewBoolVar('pair_%i_%i_%i_%i'%(x1,x2,y1,y2))
        
    for x1,x2,y1,y2,pen in d['pair']:
        model.Add(pair[(x1,x2,y1,y2)]==True).OnlyEnforceIf(same[(x1,x2)]).OnlyEnforceIf(same[(y1,y2)]).OnlyEnforceIf(same[(x1,y1)].Not())
        
    for x1,x2,y1,y2,pen in d['pair']:
        objective_terms.append((pair[(x1,x2,y1,y2)]) * pen * Wpen) #Sum of Raw Scores
    model.Maximize(sum(objective_terms))
    #model.AddHint(x[0],1)
    
    #Hint: worse results
    #for i in range(Nhint):
     #   model.AddHint(x[i],1)
	 
    print('#Start Solving...')
    # Solve
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = sec
    solver.parameters.num_search_workers = Ncpu
    solver.parameters.log_search_progress = True
    solution_printer = cp_model.ObjectiveSolutionPrinter()
    #status = solver.Solve(model)
    status = solver.SolveWithSolutionCallback(model, solution_printer)
    results=[]
    if status == pywraplp.Solver.OPTIMAL or pywraplp.Solver.FEASIBLE:
        for f in range(max_fid+1):
            for c in range(max_cid+1):
                if solver.BooleanValue(assign[(f,c)]):
                    print(f,c,'ID=',idtbl[(f,c)]) #order->ID
                    results.append(d['fragment'][idtbl[(f,c)]])
    #print(results)
    return results, solver.ObjectiveValue()


# In[7]:


def ShowPDB(results,d,filename):
    chainIdx=string.ascii_uppercase + string.ascii_lowercase + string.ascii_uppercase + string.ascii_lowercase #Max 24*2
    outlines=''
    #print(results)
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
            #[mid,fid,cid,rawsco,zsco,ali]
            ali = [int(i) for i in frg[5].split(',')[0:-1]]
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


# In[8]:


#=========================================================
def main():
    args = get_args()
    """Entry point of the program."""
    # Instantiate the data problem.
    filename = args.LIST
    #filename = 'TMP.txt'
    #filename = 'MTX_vrp_l5.txt'
    data = create_data_model(filename)
    print('#Nfragments=',len(data['fragment']))
    print('#Nnode=',len(data['node']))

    #Assemble Part
    Nsub = args.Nstock
    Ncpu = args.Ncpu
    Sec = args.SecCp
    Niter = args.Niter
    #Weight
    Weight = args.Weight
    Fkeep = args.Keep

    #Update penalty
    for i in range(len(data['pair'])):
        data['pair'][i][4] = int(data['pair'][i][4] * Weight)

    TMP_FRAG = data['fragment']
    TMP_SOLUTION={}
    TMP_SOLUTION['results']=[]
    TMP_SOLUTION['score']=0
    
    tmp_sol,score=Assemble(data,Sec,Ncpu,Fkeep)

    TMP_SOLUTION['results']=tmp_sol
    TMP_SOLUTION['score'] = score


    print("##BEstScore=",TMP_SOLUTION['score'],len(TMP_SOLUTION['results']))
    #Make Final Model
    outputfile = args.OutNamePath
    ShowPDB(TMP_SOLUTION['results'],data,outputfile)

if __name__ == '__main__':
    main()

