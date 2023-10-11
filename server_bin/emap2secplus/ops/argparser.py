#
# Copyright (C) 2020 Xiao Wang
# Email:xiaowang20140001@gmail.com wang3702@purdue.edu
#

from collections import defaultdict
import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-O',type=str, help='Output file path. (str)')
    parser.add_argument('-F',type=str, help='Input map file path. (str)')
    parser.add_argument('-M', type=str,  default="best_model",help='Pre-trained model path.  (str) Default value: "best_model"')
    parser.add_argument('--mode',type=int,required=True,help='Control Mode for program: 0: cryo_READ structure modeling. Required parameter. (Integer), Default value: 0')
    parser.add_argument("--type",type=int,required=True,help="Control type of predictions: 0: atom prediction, 1: amino acid predictions")
    parser.add_argument("--contour",type=float,default=0,help="Contour level for input map, suggested 0.5*[author_contour]. (Float), Default value: 0.0")
    parser.add_argument("--stride",type=int,default=8,help="Stride for scanning of deep learning model. (Integer), Default value: 16.")
    parser.add_argument("--box_size",type=int,default=32,help="Input box size for deep learning model. (Integer), Default value: 64")
    parser.add_argument("--gpu",type=str,default=None,help="Specify the gpu we will use. (str), Default value: None.")
    parser.add_argument('--batch_size', type=int, default=16, help='Batch size for inference of network. (Integer), Default value: 8.')
    parser.add_argument("--num_workers",type=int,default=4,help="number of workers to fetch data for GPU inference. (Integer), Default value: 4")
    args = parser.parse_args()
    params = vars(args)
    return params
