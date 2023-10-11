
import numpy as np
import os

def gen_input_data(map_data,voxel_size,stride,contour,train_save_path):
    scan_x, scan_y, scan_z = map_data.shape
    count_voxel = 0
    Coord_Voxel = []
    for x in range(0, scan_x, stride):
        x_end = min(x + voxel_size, scan_x)
        for y in range(0, scan_y, stride):
            y_end = min(y + voxel_size, scan_y)
            for z in range(0, scan_z, stride):
                z_end = min(z + voxel_size, scan_z)
                if x_end < scan_x:
                    x_start = x
                else:
                    x_start = x_end - voxel_size

                    if x_start<0:
                        x_start=0
                if y_end < scan_y:
                    y_start = y
                else:
                    y_start = y_end - voxel_size

                    if y_start<0:
                        y_start=0
                if z_end < scan_z:
                    z_start = z
                else:
                    z_start = z_end - voxel_size

                    if z_start<0:
                        z_start=0
                #already normalized
                segment_map_voxel = np.zeros([voxel_size,voxel_size,voxel_size])
                segment_map_voxel[:x_end-x_start,:y_end-y_start,:z_end-z_start]=map_data[x_start:x_end, y_start:y_end, z_start:z_end]
                if contour<=0:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel>0))
                    meaningful_density_ratio = meaningful_density_count/float(voxel_size**3)
                    if meaningful_density_ratio<=0.001:
                        print("meaningful density ratio %f, skip it!"%meaningful_density_ratio)
                        continue
                else:
                    meaningful_density_count = len(np.argwhere(segment_map_voxel > contour))
                    meaningful_density_ratio = meaningful_density_count / float(voxel_size ** 3)
                    if meaningful_density_ratio <= 0.001:
                        print("no meaningful density ratio, skip it!")
                        continue
                cur_path = os.path.join(train_save_path,"input_"+str(count_voxel)+".npy")
                np.save(cur_path,segment_map_voxel)
                Coord_Voxel.append([x_start,y_start,z_start])
                count_voxel+=1
    Coord_Voxel = np.array(Coord_Voxel)
    coord_path = os.path.join(train_save_path,"Coord.npy")
    np.save(coord_path,Coord_Voxel)
    return Coord_Voxel
