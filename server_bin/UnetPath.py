import os
import subprocess
import numpy as np
import Bio
import csv
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

MaxDis = 38820

def get_args(): 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("NODE", type=str, help="NODE file (in PDB format)")
    parser.add_argument('--Nchain', type=int, dest='Nchain', default=10, help='Number of Chains')
    parser.add_argument('--SecTrace', type=int, dest='SecTrace', default=300, help='Computatinal Time Limit for Tracing')
    parser.add_argument('--SecAssemble', type=int, dest='SecCp', default=300, help='Computatinal Time Limit for Assemble')
    parser.add_argument('--Nalignment', type=int, dest='Nali', default=20, help='Number of alignments for each path')
    #parser.add_argument('--OutThread', type=str, dest='OutNameThread', default='thread.pdb', help='Final Outpt PDB file')
    parser.add_argument('--OutCSV', type=str, dest='OutNamePath', default='path.csv', help='Final Outpt PDB file')
    args = parser.parse_args()
    return args

def create_data_model(filename):
    """Stores the data for the problem."""
    data = {}
    node=[]
    dmtx=[]
    smtx=[]
    seq=[]
    with open(filename) as f:
        for l in f:
            if l.startswith('DMTX'):
                x = int(l.split()[1])
                y = int(l.split()[2])
                d = int(l.split()[3])
                min_prob = int(l.split()[4]) #lowest Prob,
                ave_prob = int(l.split()[5])
                dmtx.append([x,y,d,min_prob,ave_prob])
            
            if l.startswith('SCO'):
                x = int(l.split()[1]) #seq
                y = int(l.split()[2]) #node
                atom_sco = float(l.split()[3]) #AAscore
                aa_sco = float(l.split()[4]) #AAscore
                smtx.append([x,y,atom_sco,aa_sco])
            
            if l.startswith('SEQ'):
                x = int(l.split()[1])
                seq.append([x,l.split()[3],int(l.split()[2])])
                
            if l.startswith('NODE'):
                x = int(l.split()[1])
                node.append([x,float(l.split()[2]),float(l.split()[3]),float(l.split()[4])])
    #Convert
    Ntot = node[-1][0] #Number of node
    Lseq = seq[-1][0] #Length of  sequence
    print('Ntot=',Ntot)
    print('Lseq=',Lseq)
    data['sequence'] = [ ['XXX',-1] for i in range(Lseq+1)]
    data['node'] = [ [0.0,0.0,0.0] for x in range(Ntot+1)]
    for row in seq:
        data['sequence'][row[0]]=[row[1],row[2]]
    for row in node:
        data['node'][row[0]]=[row[1],row[2],row[3]]
    
    #Score
    data['score_matrix'] = [[ 0 for y in range(Ntot+1)] for x in range(Lseq+1)]
    data['score_matrix_nonegative'] = [[ 0 for y in range(Ntot+1)] for x in range(Lseq+1)]
    for row in smtx:
        data['score_matrix'][row[0]][row[1]]=int(row[3]*100)
        data['score_matrix_nonegative'][row[0]][row[1]]=int(row[3]*100)
        if data['score_matrix_nonegative'][row[0]][row[1]] < 1:
            data['score_matrix_nonegative'][row[0]][row[1]] = 1

    #initial data
    data['dmtx'] = [[ [0,0,0] for x in range(Ntot+1)] for y in range(Ntot+1)]

    for x,y,d,min_prob,ave_prob in dmtx:
        data['dmtx'][x][y]=[d,min_prob,ave_prob]

    data['num_vehicles'] = 1
    data['depot'] = 0
    return data

def print_solution(data, manager, routing, solution):
    """Prints solution on console."""
    print(f'Objective: {solution.ObjectiveValue()}')
    max_route_distance = 0
    all_path = []
    for vehicle_id in range(data['num_vehicles']):
        tmp_path = []
        index = routing.Start(vehicle_id)
        plan_output = 'Route for vehicle {}:\n'.format(vehicle_id)
        route_distance = 0
        while not routing.IsEnd(index):
            plan_output += ' {} -> '.format(manager.IndexToNode(index))
            previous_index = index
            index = solution.Value(routing.NextVar(index))
            route_distance += routing.GetArcCostForVehicle(
                previous_index, index, vehicle_id)
            if manager.IndexToNode(index) != 0:
                tmp_path.append(manager.IndexToNode(index))
        plan_output += '{}\n'.format(manager.IndexToNode(index))
        plan_output += 'Distance of the route: {}m\n'.format(route_distance)
        print(plan_output)
        max_route_distance = max(route_distance, max_route_distance)
        all_path.append(tmp_path)
    print('Maximum of the route distances: {}A'.format(max_route_distance))
    return all_path

def print_solution_csv_0(data,p,outname):
    chainIdx=string.ascii_uppercase + string.ascii_lowercase + string.ascii_uppercase + string.ascii_lowercase #Max 24*2

    """Prints solution on console."""
    Natm=1
    plan_output = ''
    ND_ARRAY=[]
    PATH_ARRAY=[]
    for AtmIdx in range(len(p)):
        #print(AtmIdx,p[AtmIdx][0],p[AtmIdx][1],p[AtmIdx][2])
        ND_ARRAY.append(['NODE',AtmIdx,p[AtmIdx][0],p[AtmIdx][1],p[AtmIdx][2]])

    for vehicle_id in range(len(data)):
        AtmList=[]
        
        #plan_output +='PATH {:7d}'.format(vehicle_id)
        for AtmIdx in data[vehicle_id]:
            #plan_output +='{:7d}'.format(AtmIdx)
            #AtmIdx = str(AtmIdx)
            #plan_output +='ATOM{:7d}  CA  {:3} {:1}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}\n'.format(Natm,'ALA',chainIdx[vehicle_id],int(Natm),p[AtmIdx][0],p[AtmIdx][1],p[AtmIdx][2],1.0)
            #Natm = Natm+1
            AtmList.append(AtmIdx)
        #plan_output += '\n'
        PATH_ARRAY.append(['PATH',vehicle_id]+AtmList)
        
    #outname = 'path.pdb'
    #Save as CSV
    with open(outname,'w') as f:
        writer = csv.writer(f)
        writer.writerows(ND_ARRAY)
        writer.writerows(PATH_ARRAY)
    #with open(outname,'w') as out:
    #    out.write(plan_output)
        
def print_solution_csv(data,p,smtx,outname):
    """Prints solution on console."""
    Natm=1
    plan_output = ''
    ND_ARRAY=[]
    PATH_ARRAY=[]
    SCO_ARRAY=[]
    for AtmIdx in range(len(p)):
        ND_ARRAY.append(['NODE',AtmIdx,p[AtmIdx][0],p[AtmIdx][1],p[AtmIdx][2]])

    for vehicle_id in range(len(data)):
        AtmList=[]
        
        for AtmIdx in data[vehicle_id]:
            AtmList.append(AtmIdx)
        PATH_ARRAY.append(['PATH',vehicle_id]+AtmList)
    
    for spos in range(0,len(smtx)):
        tmp=[]
        for npos in range(0,len(smtx[0])):
            tmp.append(smtx[spos][npos])
        SCO_ARRAY.append(['SCO',spos]+tmp)
    #Save as CSV
    with open(outname,'w') as f:
        writer = csv.writer(f)
        writer.writerows(ND_ARRAY)
        writer.writerows(PATH_ARRAY)
        writer.writerows(SCO_ARRAY)       


def ConvertDmtx(data,rnd_range,sig):
    #convert dmtx
    Npos = len(data['dmtx'])
    print('Npos=',Npos)
    data['distance_matrix'] = [[ 0 for x in range(Npos)] for y in range(Npos)]
    data['real_d_matrix'] = [[ 0 for x in range(Npos)] for y in range(Npos)]
    def EdgeDistance (dist,mu,sig,minP,aveP):
        if minP==0:
            return 100+dist
        disP=np.exp(-np.power(dist - mu, 2.)/(2*np.power(sig,2.)))
        #return int(100.00-(disP*aveP)) #bad ??
        return int(100.00-(disP*minP))
        
    for x in range(1,Npos):
        for y in range(x+1,Npos):
            #[d,min_prob,ave_prob]
            d,min_prob,ave_prob = data['dmtx'][x][y]
            #Key param sig
            #sig = 10.0
            if rnd_range <= 0:
                AddRandom = 0
            else:
                AddRandom = random.randrange(2*rnd_range)-rnd_range
            if AddRandom < 0:
                AddRandom = 0
            data['distance_matrix'][x][y]=EdgeDistance(d,38,sig,min_prob,ave_prob) + AddRandom
            data['distance_matrix'][y][x]=data['distance_matrix'][x][y]
            #if min_prob > 0:
            #print(x,y,d,min_prob,ave_prob,data['distance_matrix'][x][y],EdgeDistance(d,38,sig,min_prob,ave_prob))
                
            data['real_d_matrix'][x][y]=d
            data['real_d_matrix'][y][x]=d
            
def VRP_solver(data):
    # Create the routing index manager.
    manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']),
                                           data['num_vehicles'], data['depot'])
    # Create Routing Model.
    routing = pywrapcp.RoutingModel(manager)

    # Create and register a transit callback.
    def distance_callback(from_index, to_index):
        """Returns the distance between the two nodes."""
        # Convert from routing variable Index to distance matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return data['distance_matrix'][from_node][to_node]

    transit_callback_index = routing.RegisterTransitCallback(distance_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # Add Distance constraint.
    dimension_name = 'Distance'
    routing.AddDimension(
        transit_callback_index,
        0,  # no slack
        100*10000,# vehicle maximum travel distance
        True,  # start cumul to zero
        dimension_name)
    distance_dimension = routing.GetDimensionOrDie(dimension_name)
    #distance_dimension.SetGlobalSpanCostCoefficient(100)

    #Drop penalty
    penalty = 200
    for node in range(1, len(data['distance_matrix'])):
        routing.AddDisjunction([manager.NodeToIndex(node)], penalty)

    search_parameters = pywrapcp.DefaultRoutingSearchParameters()

    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.AUTOMATIC)
    search_parameters.local_search_metaheuristic = (
        routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH)
    search_parameters.time_limit.seconds = data['CPUtime']
    search_parameters.log_search = True

    if len(data['initial_routes'])>0:
        #data['initial_routes'] = ALL_PATH
        initial_solution = routing.ReadAssignmentFromRoutes(data['initial_routes'],True)
        solution = routing.SolveFromAssignmentWithParameters(initial_solution,search_parameters)
    else:
        # Solve the problem.
        solution = routing.SolveWithParameters(search_parameters)
    #solution = routing.SolveWithParameters(search_parameters)
    # Print solution on console.
    RES_PATH=[]
    if solution:
        RES_PATH=print_solution(data, manager, routing, solution)
        #check coverage
        Npath = 0
        for path in RES_PATH:
            Npath = Npath + len(path)
        RatePath = float(Npath/(len(data['distance_matrix'])-1.00))
        print(Npath,'Rate:',RatePath)
        if RatePath < 0.9:
            print('Too Short!!')
        print('PASS')
    else:
        print('No solution found !')
    #print_solution_pdb(ALL_PATH,data['node'],'path_sig'+str(path_round)+'.pdb')
    
    return RES_PATH


def ClashCheck(data,dmtx,sdata):
    N=len(data)
    out = [[ 0 for y in range(N)] for x in range(N)]
    for i in range(N):
        seq1=[ali[0] for ali in data[i]['alignment']]
        pos1=[ali[1] for ali in data[i]['alignment']]
        
        for j in range(i+1,N):
            seq2=[ali[0] for ali in data[j]['alignment']]
            pos2=[ali[1] for ali in data[j]['alignment']]
            
            #Seq overlap
            seq_and = set(seq1) & set(seq2)
            #print('SEQAND',seq_and)
            if len(seq_and) > 1:
                out[i][j]=1
                out[j][i]=1
                continue
            seq_and_list=list(seq_and)
            if (len(seq_and)==1 and seq_and_list[0]!=-1):
                out[i][j]=1
                out[j][i]=1
                continue
            
            pos_and = set(pos1) & set(pos2)

            if len(pos_and) > 1:
                out[i][j]=1
                out[j][i]=1
                continue
            pos_and_list=list(pos_and)
            if (len(pos_and)==1 and pos_and_list[0]!=-1):
                out[i][j]=1
                out[j][i]=1
                continue
            
            #Connection check
            #chain-id
            if sdata[seq1[0]][1] == sdata[seq2[0]][1]:
                Diff1_2 = abs(seq1[-1] - seq2[0])
                if dmtx[pos1[-1]][pos2[0]] > 50 * Diff1_2:
                    out[i][j]=1
                    out[j][i]=1
                    #print(i,j,Diff1_2,dmtx[pos1[-1]][pos2[0]])
                    continue

                Diff2_1 = abs(seq2[-1] - seq1[0])
                if dmtx[pos1[0]][pos2[-1]] > 50 * Diff2_1:
                    out[i][j]=1
                    out[j][i]=1
                    #print(i,j,Diff2_1,dmtx[pos1[0]][pos2[-1]])
                    continue
                
    return out

def SimpleDP(p,dmtx,smtx,refmtx,mode,Niter,Kca):
    RESULTs=[]
    GapNodePen = -1000 #missing Node
    GapAAPen = -1000 #missing AA
    #Clash_Pen = -200 #needs refine!
    pth = p
    if mode == 1:
        pth = [i for i in reversed(p)]
    print('Path',pth)
    Npos = len(pth)+1
    Naa  = len(smtx)
    SMTX = [[ 0 for y in range(Npos)] for x in range(Naa)]

    #Fill SMTX
    for x in range(1,Naa):
        for y in range(1,Npos):
            SMTX[x][y] = smtx[x][pth[y-1]]


    for rnd in range(Niter):
        MTX = [[ {'Dia':0,'Sco':0,'Pre':[0,0]} for y in range(Npos)] for x in range(Naa)]

        #Fill DiaMTX
        for x in range(1,Naa):
            for y in range(1,Npos):
                
                UpSco   = MTX[x][y-1]['Sco'] + GapNodePen #missing Node
                LeftSco = MTX[x-1][y]['Sco'] + GapAAPen #missing AA
                DiaSco  = MTX[x-1][y-1]['Sco'] + SMTX[x][y]
                #Check Previous aligned position
                NowPos=pth[y-1]
                PrePos=pth[MTX[x-1][y-1]['Pre'][1]-1]
                #if MTX[x-1][y-1]['Pre'][1]!=0:
                #    dist = dmtx[NowPos][PrePos]
                #    if dist < 20 or dist > 50:
                #        DiaSco = DiaSco + (dist-38)*(dist-38)*Kca

                if DiaSco <= 0 and LeftSco <= 0 and UpSco<= 0:
                    MTX[x][y]['Dia'] = 0
                elif DiaSco >=UpSco and DiaSco>=LeftSco:
                    MTX[x][y]['Sco'] = DiaSco
                    MTX[x][y]['Dia'] = 1
                    MTX[x][y]['Pre'] = [x,y]
                elif UpSco >= DiaSco and UpSco >=LeftSco:
                    MTX[x][y]['Sco'] = UpSco
                    MTX[x][y]['Dia'] = 2
                    MTX[x][y]['Pre'] = MTX[x][y-1]['Pre']
                elif LeftSco >= DiaSco and LeftSco >=UpSco:
                    MTX[x][y]['Sco'] = LeftSco
                    MTX[x][y]['Dia'] = 3
                    MTX[x][y]['Pre'] = MTX[x-1][y]['Pre']
                #print(x,y,MTX[x][y],SMTX[x][y],DiaSco,PrePos)

        #find highest
        Hsco = 0
        Hpos = []
        for x in range(1,Naa):
            for y in range(1,Npos):
                #if TmpMTX[x][y] > Hsco:
                if MTX[x][y]['Dia'] != 1:
                    continue
                if MTX[x][y]['Sco'] >= Hsco:
                    Hsco = MTX[x][y]['Sco']
                    Hpos = [x,y]

        if Hsco == 0: #No new alignments
            break
        print('TraceBack',rnd,Hsco, Hpos)
        #TraceBack....
        ALI=[]
        while True:
            #print(Hpos,MTX[Hpos[0]][Hpos[1]])
            if MTX[Hpos[0]][Hpos[1]]['Dia'] == 1:#Dia
                ALI.append([Hpos[0],Hpos[1]])
                Hpos[0] = Hpos[0] -1
                Hpos[1] = Hpos[1] -1
            elif MTX[Hpos[0]][Hpos[1]]['Dia'] == 2:#Up
                ALI.append([-1,Hpos[1]])
                Hpos[0] = Hpos[0]
                Hpos[1] = Hpos[1] -1
            elif MTX[Hpos[0]][Hpos[1]]['Dia'] == 3:#Left
                ALI.append([Hpos[0],-1])
                Hpos[0] = Hpos[0] -1
                Hpos[1] = Hpos[1]
            elif MTX[Hpos[0]][Hpos[1]]['Dia'] == 0:#Stop
                break

        #Convert...Sequence Position and Node
        OUT = []
        RawSco = 0.0
        for i in reversed(ALI):
            if i[0] > 0 and i[1] > 0:
                OUT.append([i[0],pth[i[1]-1]])
                RawSco = RawSco + refmtx[i[0]][pth[i[1]-1]]
            #Update SMTX Masking
            if i[0]>0 and i[1] > 0:
                if SMTX[i[0]][i[1]] > 0:
                    SMTX[i[0]][i[1]]=0
        print('RawSco=',RawSco)
        #print(OUT)
        if len(OUT)>3:#at least 4 residues
            #RESULTs.append({'alignment':OUT,'score':Hsco})
            RESULTs.append({'alignment':OUT,'score':RawSco})
            #print(OUT)

    return RESULTs

def ShowPDB_fragments(res,d,outname):
    chainIdx=string.ascii_uppercase + string.ascii_lowercase + string.ascii_uppercase + string.ascii_lowercase #Max 24*2
    outlines = ''
    pairs=[]
    Natm = 1
    p=d['node']
    for f in res:
        #print(f['alignment'])
        for pair in f['alignment']:
            pairs.append(pair)
    pairs.sort(key=lambda row:(row[0]))
    #print(pairs)
    
    #Fill Missing Residues
    COORD=[]
    for pair in pairs:
        #AtmIdx = str(pair[1])
        AtmIdx = pair[1]
        SeqIdx = pair[0]
        COORD.append([SeqIdx,p[AtmIdx]])
    
    LOOP=[]
    if False:
        for i in range(len(COORD)-1):
            if COORD[i][0]+1 != COORD[i+1][0]: #Missing Residues
                #print('Fill',COORD[i][0], '->',COORD[i+1][0])
                if COORD[i+1][0] - COORD[i][0] < 8: #Small Gap
                    SeqDiff = COORD[i+1][0] - COORD[i][0]
                    VEC=[0.0,0.0,0.0]
                    print(COORD[i+1][1],COORD[i][1])
                    VEC[0]=(COORD[i+1][1][0]-COORD[i][1][0])/float(SeqDiff)
                    VEC[1]=(COORD[i+1][1][1]-COORD[i][1][1])/float(SeqDiff)
                    VEC[2]=(COORD[i+1][1][2]-COORD[i][1][2])/float(SeqDiff)
                
                    for j in range(1,SeqDiff):
                        #print(COORD[i][0],j,COORD[i][1][0]+VEC[0]*j)
                        LOOP.append([COORD[i][0]+j,[COORD[i][1][0]+VEC[0]*j,COORD[i][1][1]+VEC[1]*j,COORD[i][1][2]+VEC[2]*j]])
    
    OUT=COORD+LOOP
    OUT.sort(key=lambda row:(row[0]))
    for pair in OUT:
        SeqIdx = pair[0]
        coords = pair[1]
        #print(pair)
        cid=d['sequence'][SeqIdx][1]
        #outlines +='ATOM{:7d}  CA  {:3} {:1}{:4d}'.format(Natm,d['sequence'][SeqIdx][0],'A',SeqIdx)
        outlines +='ATOM{:7d}  CA  {:3} {:1}{:4d}'.format(Natm,d['sequence'][SeqIdx][0],chainIdx[cid],SeqIdx)
        
        outlines +='    {:8.3f}{:8.3f}{:8.3f}  1.00{:6.2f}\n'.format(coords[0],coords[1],coords[2],1.0)
        Natm = Natm+1

    #outname = 'thread.pdb'
    with open(outname,'w') as out:
        out.write(outlines)



#=========================================================
def main():
    args = get_args()
    
    # Instantiate the data problem.
    filename = args.NODE
    
    data = create_data_model(filename)


    ALL_PATH=[]
    data['initial_routes']=[]
    for path_round in range(5,31):#sigma 0.5-3.0
    
        sig = float(path_round)
        ConvertDmtx(data,0,sig) #rndom range, sigma
    
        data['num_vehicles'] = args.Nchain
        data['CPUtime'] = args.SecTrace
    
        tmp_path=VRP_solver(data)
        data['initial_routes']=tmp_path
        #print_solution_pdb(tmp_path,data['node'],'path_sig'+str(path_round)+'.pdb')
        
        #chek NR
        for p1 in tmp_path:
            flag =  False
            for p2 in ALL_PATH:
                if p1 == p2:
                    print('SAME',p1,p2)
                    flag = True
                    break
                p1rv = [i for i in reversed(p1)]
                if p1rv == p2:
                    print('SAME',p1rv,p2)
                    flag = True
                    break
            if flag == False:
                ALL_PATH.append(p1)
                
    print_solution_csv(ALL_PATH,data['node'],data['score_matrix'],args.OutNamePath)

if __name__ == '__main__':
    main()
