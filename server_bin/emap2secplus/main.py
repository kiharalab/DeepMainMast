import os
from ops.argparser import argparser
from ops.os_operation import mkdir
import time
import subprocess
import time
import mrcfile
import numpy as np
def get_gpu_memory_usage():
    result = subprocess.run(['nvidia-smi', '--query-gpu=memory.used', '--format=csv,nounits,noheader'], stdout=subprocess.PIPE)
    output = result.stdout.decode('utf-8').strip()
    memory_usage = [int(x) for x in output.split('\n')]
    return memory_usage
if __name__ == "__main__":
    params = argparser()
    if params['mode']==0:
        gpu_id = params['gpu']
        if gpu_id is not None:
            os.environ["CUDA_VISIBLE_DEVICES"] = gpu_id
        cur_map_path = os.path.abspath(params['F'])
        model_path = os.path.abspath(params['M'])
        save_path = f"results/{params['F'].split('/')[1]}/unet"
        save_path = os.path.abspath(save_path)
        map_name = params['F'].split('/')[1]
        from data_processing.Unify_Map import Unify_Map
        cur_map_path = Unify_Map(cur_map_path,os.path.join(save_path,map_name+"_unified.mrc"))
        from data_processing.Resize_Map import Resize_Map
        cur_map_path = Resize_Map(cur_map_path,os.path.join(save_path,map_name+".mrc"))
        from data_processing.map_utils import segment_map
        cur_new_map_path = os.path.join(save_path,map_name+"_segment.mrc")
        cur_map_path = segment_map(cur_map_path,cur_new_map_path,contour=0) #save the final prediction prob array space
#        gpu = get_gpu_memory_usage()
#        free_gpu =12288 - max(gpu)
#        while free_gpu<=8000:
#            print("waiting for available gpu")
#            time.sleep(60)
#            gpu = get_gpu_memory_usage()
#            free_gpu =12288 - max(gpu)
        with mrcfile.open(cur_map_path,permissive=True) as mrc:
            data=mrc.data
        if np.max(data)<=params['contour']:
            print("!!!Warning!!!Please provide contour level lower than maximum density value to run!")
            #exit()
            from data_processing.map_utils import automate_contour
            params['contour']=automate_contour(data)
            print("!!!Warning!!!reset contour level: %f"%params['contour'])
        from predict.detect_map import detect_map
        detect_map(cur_map_path,model_path,save_path,
                                  params['box_size'],params['stride'],
                                  params['batch_size'],params['contour'],params)



