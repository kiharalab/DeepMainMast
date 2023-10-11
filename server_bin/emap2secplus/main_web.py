import os
from ops.argparser import argparser
from ops.os_operation import mkdir
import time
import subprocess
import time
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
        save_path = f"{params['O']}/unet"
        map_name = params['F'].split('/')[-1].split('.')[-2]
        from data_processing.Resize_Map import Resize_Map
        cur_map_path = Resize_Map(cur_map_path,os.path.join(save_path,map_name+".mrc"))
        from predict.detect_map import detect_map
        gpu = get_gpu_memory_usage()
        free_gpu =12288 - max(gpu)
        while free_gpu<=8000:
            print("waiting for available gpu")
            time.sleep(60)
            gpu = get_gpu_memory_usage()
            free_gpu =12288 - max(gpu)
        detect_map(cur_map_path,model_path,save_path,
                                  params['box_size'],params['stride'],
                                  params['batch_size'],params['contour'],params)



