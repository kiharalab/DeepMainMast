
import os
import numpy as np
from data_processing.gen_input_data import gen_input_data
from data_processing.Single_Dataset import Single_Dataset
import torch
import torch.nn as nn
from model.Small_Unet_3Plus_DeepSup import Small_UNet_3Plus_DeepSup
from predict.make_predictions import make_predictions

def unet_detect_protein(map_data,resume_model_path,voxel_size,
                    stride,batch_size,train_save_path,contour,params,original_map_path,save_root):
    coord_path = os.path.join(train_save_path, "Coord.npy")
    if os.path.exists(coord_path):
        Coord_Voxel = np.load(coord_path)
    else:
        Coord_Voxel = gen_input_data(map_data, voxel_size, stride, contour, train_save_path)
    overall_shape = map_data.shape
    test_dataset = Single_Dataset(train_save_path, "input_")
    test_loader = torch.utils.data.DataLoader(
        test_dataset,
        pin_memory=True,
        batch_size=batch_size,
        shuffle=False,
        num_workers=params['num_workers'],
        drop_last=False)
    if params['type']==0:
        base_class = 7
        refer_name = "atom"
        label_list=['BG','N',"CA","C","O","CB","Others"]
        pre_name="atom"
    elif params['type']==1:
        base_class = 20
        refer_name = "sigmoidAA"
        #check the residue types
        label_list = ["ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU", "ASP", "GLU", "LYS", "ARG", "SER", "THR", "TYR",
                "HIS", "CYS", "ASN", "TRP", "GLN", "GLY"]
        pre_name = "sigmoidAA"
    else:
        print("only support --type 0 or 1. %d type is not suppored"%params['type'])
        exit()

    model = Small_UNet_3Plus_DeepSup(in_channels=1,
                                     n_classes=base_class,
                                     feature_scale=4,
                                     is_deconv=True,
                                     is_batchnorm=True)
    model = model.cuda()
    model = nn.DataParallel(model, device_ids=None)
    state_dict = torch.load(resume_model_path)
    msg = model.load_state_dict(state_dict['state_dict'])
    print("model loading: ", msg)
    #cur_prob_path = os.path.join(train_save_path, refer_name + "_predictprob.npy")
    # #cur_label_path = os.path.join(train_save_path, refer_name + "_predict.npy")
    # if os.path.exists(cur_prob_path) and os.path.exists(cur_label_path):
    #     Prediction_Matrix = np.load(cur_prob_path)
    #     Prediction_Label = np.load(cur_label_path)
    #
    # else:
    make_predictions(test_loader, model, Coord_Voxel,
                        voxel_size, overall_shape,
                         base_class,save_root,
                     pre_name,label_list,original_map_path,run_type=params['type'])

        # np.save(cur_prob_path, Prediction_Matrix)
        # np.save(cur_label_path, Prediction_Label)
    # return Prediction_Matrix
