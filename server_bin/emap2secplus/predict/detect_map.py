


from ops.os_operation import mkdir

import os
import mrcfile
import numpy as np
from data_processing.map_utils import save_predict_specific_map,permute_ns_coord_to_pdb,find_top_density
from predict.unet_detect_protein import unet_detect_protein

def detect_map(input_map_path,resume_model_path,save_path,
                              voxel_size,stride,batch_size,contour,params):
    """
    :param input_map_path:
    :param resume_model_path:
    :param voxel_size:
    :param stride:
    :param batch_size:
    :return:
    """

    print("TRYING TO OPEN MAP")
    with mrcfile.open(input_map_path, permissive=True) as map_mrc:
         #normalize data
        map_data = np.array(map_mrc.data)
        # get the value serve as 1 in normalization
        map_data[map_data < 0] = 0
        print("map density range: %f %f"%(0,np.max(map_data)))
        percentile_98 = find_top_density(map_data,0.98)

        print("HAHA")
        print("map hist log percentage 98: ",percentile_98)
        map_data[map_data > percentile_98] = percentile_98
        min_value = np.min(map_data)
        max_value = np.max(map_data)
        map_data = (map_data-min_value)/(max_value-min_value)
        nxstart, nystart, nzstart = map_mrc.header.nxstart, \
                                    map_mrc.header.nystart, \
                                    map_mrc.header.nzstart
        orig = map_mrc.header.origin
        orig = str(orig)
        orig = orig.replace("(", "")
        orig = orig.replace(")", "")
        orig = orig.split(",")
        nstart = [nxstart, nystart, nzstart]
        mapc = map_mrc.header.mapc
        mapr = map_mrc.header.mapr
        maps = map_mrc.header.maps
        print("detected mode mapc %d, mapr %d, maps %d" % (mapc, mapr, maps))
        nstart = permute_ns_coord_to_pdb(nstart, mapc, mapr, maps)
        new_origin = []
        for k in range(3):
            new_origin.append(float(orig[k]) + float(nstart[k]))

        print("Origin:", new_origin)
        train_save_path = os.path.join(save_path,"Input")
        mkdir(train_save_path)
        #adjust the contour level by the maximum value
        print("given contour %f"%contour)
        contour = contour/percentile_98
        print("revised contour %f"%contour)
        detection_protein = unet_detect_protein(map_data, resume_model_path,
                                           voxel_size, stride, batch_size,
                                           train_save_path, contour,
                                    params,input_map_path,save_path)

        # detection_protein = unet_detect_protein(map_data, resume_model_path,
        #                                   voxel_size, stride, batch_size,
        #                                   train_save_path, contour, params)


        if detection_protein is not None:
            if params['type']==0:
                label_list=['BG','N',"CA","C","O","CB","Others"]
                pre_name="atom"
            elif params['type']==1:
                #check the residue types
                label_list = ["ALA", "VAL", "PHE", "PRO", "MET", "ILE", "LEU", "ASP", "GLU", "LYS", "ARG", "SER", "THR", "TYR",
                        "HIS", "CYS", "ASN", "TRP", "GLN", "GLY"]
                pre_name = "sigmoidAA"
            for k, base_name in enumerate(label_list):
                cur_map_path = os.path.join(save_path, pre_name+"_" + str(base_name) + ".mrc")
                save_predict_specific_map(cur_map_path, k, detection_protein, input_map_path)

    print("finish prediction")
