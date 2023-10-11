import os
import sys
import csv

job_dir = sys.argv[1]
out_dir = sys.argv[2]

def get_daq(fname):
    ave_aa=0.0
    ave_ca=0.0
    tot_aa=0.0
    tot_ca=0.0
    with open(fname) as f:
        for l in f:
            if l.startswith('AVG   Calpha'):
                ave_ca=float(l.split()[4])
                ave_aa=float(l.split()[6])
            if l.startswith('TOTAL Calpha'):
                tot_ca=float(l.split()[4])
                tot_aa=float(l.split()[6])

    return ave_ca,tot_ca,ave_aa,tot_aa

def get_dot(fname):
    vol=0.0
    dot=0.0
    with open(fname) as f:
        for l in f:
            if l.startswith('Overlap'):
                vol=float(l.split()[2].split('/')[1])
                dot=float(l.split()[10])
    if vol==0.0:
        return 0,0
    return dot,dot/vol

methods = ["all", "VESPER", "DMonly", "AFonly"]
models_dict = {}

#tbl=[['model','DAQ','DOT','TOTAL']]
tbl=[]
for method in methods:
    for mid in [1,2,3,4,5]: #5 models
        model_path = f"{job_dir}/CM_{method}/S_singletgt_000{mid}."
        if os.path.exists(model_path + "daq"):
            daq_score = get_daq(model_path + "daq")[2]  
            dot_score = get_dot(model_path + "dot")[1]  
            tbl.append([model_path+"pdb",daq_score+dot_score,daq_score,dot_score])

models = sorted(tbl, key=lambda x: x[1], reverse=True)

#Save csv

rank = 1
for model_path, score,daq,dot in models:
    print(model_path, score)
    daq_file=model_path.replace('.pdb','.daq')
    os.system(f"cp {model_path} {out_dir}/rank{rank}.pdb")
    os.system(f"cp {daq_file} {out_dir}/rank{rank}.daq")
    rank += 1

csvfile=f"{out_dir}/ranking.csv"
models = [['MODEL','TOTAL','DAQ','DOT']] + models
with open(csvfile,mode="w") as file:
    writer=csv.writer(file)
    writer.writerows(models)

