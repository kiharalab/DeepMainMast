
import numpy as np
import os
import datetime
import time
import torch
import torch.nn as nn
from ops.Logger import AverageMeter,ProgressMeter
from ops.os_operation import mkdir
from data_processing.map_utils import save_predict_specific_map
def make_predictions_old(test_loader,model,Coord_Voxel,voxel_size,overall_shape,
                     num_classes,run_type=0):
    avg_meters = {'data_time': AverageMeter('data_time'),
                  'train_time': AverageMeter('train_time')}
    progress = ProgressMeter(
        len(test_loader),
        [avg_meters['data_time'],
         avg_meters['train_time'],
         ],
        prefix="#Eval:")
    model.eval()
    end_time = time.time()
    scan_x, scan_y, scan_z = overall_shape
    Prediction_Matrix = np.zeros([num_classes,overall_shape[0],overall_shape[1],overall_shape[2]])
    Count_Matrix = np.zeros(overall_shape)
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            # input, atom_target, nuc_target = data
            input, cur_index = data
            #print(input.shape)
            avg_meters['data_time'].update(time.time() - end_time, input.size(0))
            cur_id = cur_index.detach().cpu().numpy()#test_loader.dataset.id_list[cur_index.detach().numpy()]
            input = input.cuda()
            outputs = model(input)
            if run_type==2:
                final_output = torch.softmax(torch.sigmoid(outputs[0]),dim=1).detach().cpu().numpy()
            elif run_type==1:
                final_output = torch.sigmoid(outputs[0]).detach().cpu().numpy()
            else:
                final_output = torch.softmax(outputs[0],dim=1).detach().cpu().numpy()
            avg_meters['train_time'].update(time.time() - end_time, input.size(0))
            progress.display(batch_idx)
            for k in range(len(cur_id)):
                tmp_index = cur_id[k]
                x_start, y_start, z_start = Coord_Voxel[int(tmp_index)]
                x_end,y_end,z_end = x_start+voxel_size,y_start+voxel_size,z_start+voxel_size
                if x_end < scan_x:
                    x_start = x_start
                else:
                    x_end = scan_x
                    x_start = x_end - voxel_size
                    if x_start<0:
                        x_start=0
                if y_end < scan_y:
                    y_start = y_start
                else:
                    y_end = scan_y
                    y_start = y_end - voxel_size
                    if y_start<0:
                        y_start=0
                if z_end < scan_z:
                    z_start = z_start
                else:
                    z_end=scan_z
                    z_start = z_end - voxel_size
                    if z_start<0:
                        z_start=0
                #print(final_output[k].shape)
                #print(Prediction_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end].shape)
                pred_label=np.argmax(final_output[k],axis=1)
                count_positive= len(np.argwhere(pred_label!=0))
                print("%d example with %d positive predictions"%(k,count_positive))
                Prediction_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end] += final_output[k][:,:x_end-x_start,:y_end-y_start,:z_end-z_start]

                Count_Matrix[x_start:x_end,y_start:y_end,z_start:z_end]+=1
            if batch_idx%1000==0:
                for j in range(num_classes):
                    count_positive = len(np.argwhere(Prediction_Matrix[j]>=0.5))
                    print("%d classes already detected %d voxels"%(j,count_positive))
            end_time = time.time()
    Prediction_Matrix = Prediction_Matrix/Count_Matrix
    #replace nan with 0
    Prediction_Matrix[np.isnan(Prediction_Matrix)] = 0
    Prediction_Label = np.argmax(Prediction_Matrix,axis=0)
    return Prediction_Matrix,Prediction_Label


def make_predictions(test_loader,model,Coord_Voxel,voxel_size,overall_shape,
                     num_classes,save_root,pre_name,label_list,original_map_path,run_type=0):
    #to be more friendly to memory
    avg_meters = {'data_time': AverageMeter('data_time'),
                  'train_time': AverageMeter('train_time')}
    progress = ProgressMeter(
        len(test_loader),
        [avg_meters['data_time'],
         avg_meters['train_time'],
         ],
        prefix="#Eval:")
    model.eval()
    end_time = time.time()
    scan_x, scan_y, scan_z = overall_shape
    output_dir = os.path.join(save_root,"Output")
    mkdir(output_dir)
    coord_map_dict={}
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            # input, atom_target, nuc_target = data
            input, cur_index = data
            #print(input.shape)
            avg_meters['data_time'].update(time.time() - end_time, input.size(0))
            cur_id = cur_index.detach().cpu().numpy()#test_loader.dataset.id_list[cur_index.detach().numpy()]
            input = input.cuda()
            outputs = model(input)
            if run_type==2:
                final_output = torch.softmax(torch.sigmoid(outputs[0]),dim=1).detach().cpu().numpy()
            elif run_type==1:
                final_output = torch.sigmoid(outputs[0]).detach().cpu().numpy()
            else:
                final_output = torch.softmax(outputs[0],dim=1).detach().cpu().numpy()
            cur_save_path = os.path.join(output_dir,"predict_%d.npy"%batch_idx)
            np.save(cur_save_path,final_output)
            avg_meters['train_time'].update(time.time() - end_time, input.size(0))
            progress.display(batch_idx)

            coord_map_dict[cur_save_path]=cur_id
    for kk in range(num_classes):
        Prediction_Matrix = np.zeros(overall_shape)
        Count_Matrix = np.zeros(overall_shape)
        for file_path in coord_map_dict:
            cur_id = coord_map_dict[file_path]
            final_output = np.load(file_path)
            for k in range(len(cur_id)):
                tmp_index = cur_id[k]
                x_start, y_start, z_start = Coord_Voxel[int(tmp_index)]
                x_end,y_end,z_end = x_start+voxel_size,y_start+voxel_size,z_start+voxel_size
                if x_end < scan_x:
                    x_start = x_start
                else:
                    x_end = scan_x
                    x_start = x_end - voxel_size
                    if x_start<0:
                        x_start=0
                if y_end < scan_y:
                    y_start = y_start
                else:
                    y_end = scan_y
                    y_start = y_end - voxel_size
                    if y_start<0:
                        y_start=0
                if z_end < scan_z:
                    z_start = z_start
                else:
                    z_end=scan_z
                    z_start = z_end - voxel_size
                    if z_start<0:
                        z_start=0
                #print(final_output[k].shape)
                #print(Prediction_Matrix[:,x_start:x_end,y_start:y_end,z_start:z_end].shape)
                # pred_label=np.argmax(final_output[k],axis=1)
                # count_positive= len(np.argwhere(pred_label!=0))
                # print("%d example with %d positive predictions"%(k,count_positive))
                Prediction_Matrix[x_start:x_end,y_start:y_end,z_start:z_end] += final_output[k][kk,:x_end-x_start,:y_end-y_start,:z_end-z_start]

                Count_Matrix[x_start:x_end,y_start:y_end,z_start:z_end]+=1
        Prediction_Matrix = Prediction_Matrix/Count_Matrix
        #replace nan with 0
        Prediction_Matrix[np.isnan(Prediction_Matrix)] = 0
        base_name = label_list[kk]
        cur_map_path = os.path.join(save_root, pre_name+"_" + str(base_name) + ".mrc")
        tmp_prediction_matrix=[Prediction_Matrix]
        save_predict_specific_map(cur_map_path, 0,tmp_prediction_matrix, original_map_path)#use 2 since prob always <=1


