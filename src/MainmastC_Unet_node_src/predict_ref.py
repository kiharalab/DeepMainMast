import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import os
from os.path import join
from torch.autograd import Variable
from torch import FloatTensor, LongTensor
import argparse
import random
from sklearn.preprocessing import OneHotEncoder

from models.reference import NeuralNet

onehot_encoder = OneHotEncoder(sparse=False, categories = [['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']])
onehot_encoder.fit([['A'], ['R'], ['N'], ['D'], ['C'], ['Q'], ['E'], ['G'], ['H'], ['I'], ['L'], ['K'], ['M'], ['F'], ['P'], ['S'], ['T'], ['W'], ['Y'], ['V']])

def read_target_list(list_fn):
	targets = []
	f = open(list_fn, 'r')
	for line in f:
	    targets.append(line.strip())
	f.close()
	return targets

def get_seq(target,seq_dir):

    target_seq_path = join(seq_dir,target+'.fasta')
    f = open(target_seq_path,'r')
    f.readline()
    target_seq = "".join(s.strip() for s in f.readlines())
    f.close()
    return target_seq.strip()

def get_feature(mode,n,i,aa_i,j= None,aa_j = None):

	if mode == "dist":
		f_i = float(i)/600
		f_j = float(j)/600
		f_n = float(n)/600
		if aa_i == 'G' or aa_j == 'G':
			is_glycine = 1
		else:
			is_glycine = 0
		return Variable(FloatTensor([f_i,f_j,f_n,is_glycine]))

	elif mode == "ori":
		f_i = float(i)/600
		f_j = float(j)/600
		f_n = float(n)/600
		aa_i_encoding = onehot_encoder.transform([[aa_i]])[0]
		aa_j_encoding = onehot_encoder.transform([[aa_j]])[0]
		return Variable(FloatTensor(np.concatenate((np.array([f_i,f_j,f_n]), np.array(aa_i_encoding), np.array(aa_j_encoding)), axis=0)))

	elif mode == "angle":
		f_i = float(i)/600
		f_n = float(n)/600
		aa_i_encoding = onehot_encoder.transform([[aa_i]])[0]
		return Variable(FloatTensor(np.concatenate((np.array([f_i,f_n]), np.array(aa_i_encoding)), axis=0)))

def load_ref_model(model_path, input_size, out_size, dropout=0.0):
	model = NeuralNet(in_size=input_size, out_size=out_size, dropout=dropout)
	model.load_state_dict(torch.load(model_path,map_location=torch.device('cpu')))
	model.eval()
	return model

def count_to_prob(count_vec):
	s = sum(count_vec)
	prob_vec = [c/s for c in count_vec]

	return prob_vec

def symmetric_mtx(mtx):
	sym_mtx = 0.5 * (mtx + mtx.transpose(0,2,1))
	return sym_mtx

def swap(mtx):
	#From L,L,C to C,L,L
	#L,L,C
	mtx = np.swapaxes(mtx,1,2)
	#L,C,L
	mtx = np.swapaxes(mtx,0,1)
	return mtx

class Pred():

	def __init__(self):

		self.dist_bins = 20
		self.omega_bins = 25
		self.theta_bins = 25
		self.ori_phi_bins = 13
		self.phi_bins = 36
		self.psi_bins = 36
		self.hbond_bins = 38
		self.sidechain_bins = 38

	def target_data(self, target, seq_dir=None):
		self.target = target
		self.sequence = get_seq(self.target,seq_dir)
		self.len = len(self.sequence)

	def get_ref_pred(self, dist_ref_model_path, omega_ref_model_path, theta_ref_model_path, ori_phi_ref_model_path, phi_ref_model_path, psi_ref_model_path, hbond_ref_model_path, sidechain_ref_model_path):

		self.dist_ref_model = load_ref_model(dist_ref_model_path, input_size=4, out_size=20)
		self.omega_ref_model = load_ref_model(omega_ref_model_path, input_size=43, out_size=25)
		self.theta_ref_model = load_ref_model(theta_ref_model_path, input_size=43, out_size=25)
		self.ori_phi_ref_model = load_ref_model(ori_phi_ref_model_path, input_size=43, out_size=13)
		self.phi_ref_model = load_ref_model(phi_ref_model_path, input_size=22, out_size=36)
		self.psi_ref_model = load_ref_model(psi_ref_model_path, input_size=22, out_size=36)
		self.hbond_ref_model = load_ref_model(hbond_ref_model_path, input_size=4, out_size=38)
		self.sidechain_ref_model = load_ref_model(sidechain_ref_model_path, input_size=4, out_size=38)
		
		length = self.len
		dist_ref_pred = np.zeros((length, length,self.dist_bins))
		omega_ref_pred = np.zeros((length, length,self.omega_bins))
		theta_ref_pred = np.zeros((length, length, self.theta_bins))
		ori_phi_ref_pred = np.zeros((length, length,self.ori_phi_bins))
		phi_ref_pred = np.zeros((length, self.phi_bins))
		psi_ref_pred = np.zeros((length, self.psi_bins))
		hbond_ref_pred = np.zeros((length, length, self.hbond_bins))
		sidechain_ref_pred = np.zeros((length, length, self.sidechain_bins))

		for i in range(length):
			for j in range(length):

				dist_feature = get_feature("dist",n=length,i=i,aa_i=self.sequence[i],j=j,aa_j=self.sequence[j])
				dist_outs = F.softmax(self.dist_ref_model(dist_feature),dim=0)
				dist_ref_pred[i][j][:] = np.array(dist_outs.data)

				#ori_feature = get_feature("ori",n=length,i=i,aa_i=self.sequence[i],j=j,aa_j=self.sequence[j])
				#omega_outs = F.softmax(self.omega_ref_model(ori_feature),dim=0)
				#omega_ref_pred[i][j][:] = np.array(omega_outs.data)
				#theta_outs = F.softmax(self.theta_ref_model(ori_feature),dim=0)
				#theta_ref_pred[i][j][:] = np.array(theta_outs.data)
				#ori_phi_outs = F.softmax(self.ori_phi_ref_model(ori_feature),dim=0)
				#ori_phi_ref_pred[i][j][:] = np.array(ori_phi_outs.data)

				#hbond_outs = F.softmax(self.hbond_ref_model(dist_feature),dim=0)
				#hbond_ref_pred[i][j][:] = np.array(hbond_outs.data)

				#sidechain_outs = F.softmax(self.sidechain_ref_model(dist_feature),dim=0)
				#sidechain_ref_pred[i][j][:] = np.array(sidechain_outs.data)

			#angle_feature = get_feature("angle",n=length,i=i,aa_i=self.sequence[i])
			#phi_outs = F.softmax(self.phi_ref_model(angle_feature),dim=0)
			#phi_ref_pred[i][:] = np.array(phi_outs.data)
			#psi_outs = F.softmax(self.psi_ref_model(angle_feature),dim=0)
			#psi_ref_pred[i][:] = np.array(psi_outs.data)

		ref_data = {}
		ref_data["dist_ref_pred"] = symmetric_mtx(swap(dist_ref_pred))
		#ref_data["omega_ref_pred"] = symmetric_mtx(swap(omega_ref_pred))
		#ref_data["theta_ref_pred"] = swap(theta_ref_pred)
		#ref_data["ori_phi_ref_pred"] = swap(ori_phi_ref_pred)
		#ref_data["hbond_ref_pred"] = swap(hbond_ref_pred)
		#ref_data["sidechain_ref_pred"] = symmetric_mtx(swap(sidechain_ref_pred))
		#ref_data["phi_ref_pred"] = phi_ref_pred.transpose()
		#ref_data["psi_ref_pred"] = psi_ref_pred.transpose()

		return ref_data
		
	def pred_dist(self,dist_pred_dir):
		p = np.load(join(dist_pred_dir,self.target+'_prediction.npz'))
		dist_model_data = {}
		dist_model_data["dist_pred"] = p["dist"]
		dist_model_data["omega_pred"] = p["mu"]#p["omega"]
		dist_model_data["theta_pred"] = p["theta"]
		dist_model_data["orientation_phi_pred"] = p["rho"]#p["orientation_phi"]
		dist_model_data["phi_pred"] = p["phi"]
		dist_model_data["psi_pred"] = p["psi"]
		dist_model_data["sidechain_pred"] = p["sidechain_center"]
		dist_model_data["hbond_pred"] = p["hbond"]
		return dist_model_data


	def pred_sidechain_center(self,sidechain_pred_dir):
		p = np.load(join(sidechain_pred_dir,self.target+'_sce.npy'))
		p = 0.5 * (p + p.transpose(0,2,1))
		return p

	def pred_hbond(self,hbond_pred_dir):
		p = np.load(join(hbond_pred_dir,self.target+'_n_o.npy'))
		return p

def main(target,fasta_dir,out_dir,dist_pred_dir,sidechain_pred_dir,hbond_pred_dir,dist_ref_model_path, omega_ref_model_path, theta_ref_model_path, ori_phi_ref_model_path, phi_ref_model_path, psi_ref_model_path, hbond_ref_model_path, sidechain_ref_model_path):

	pred = Pred()
	pred.target_data(target=target,seq_dir=fasta_dir)
	print("Computing reference prob for",target)
	ref_data = pred.get_ref_pred(dist_ref_model_path, omega_ref_model_path, theta_ref_model_path, ori_phi_ref_model_path, phi_ref_model_path, psi_ref_model_path, hbond_ref_model_path, sidechain_ref_model_path)
	#print("Combining all predictions into final file for",target)
	#dist_model_data = pred.pred_dist(dist_pred_dir)
	#sidechain_center_pred = pred.pred_sidechain_center(sidechain_pred_dir)
	#hbond_pred = pred.pred_hbond(hbond_pred_dir)

	s = join(out_dir, target + '_final')

	np.savez(s,
		#dist=dist_model_data["dist_pred"],
		#phi=dist_model_data["phi_pred"],
		#psi=dist_model_data["psi_pred"],
		#omega=dist_model_data["omega_pred"],
		#theta=dist_model_data["theta_pred"],
		#orientation_phi=dist_model_data["orientation_phi_pred"],
		#hbond=hbond_pred,
		#sidechain_center=sidechain_center_pred,
		#hbond=dist_model_data["hbond_pred"],
		#sidechain_center=dist_model_data["sidechain_pred"],
		dist_ref_pred=ref_data["dist_ref_pred"],
		#omega_ref_pred=ref_data["omega_ref_pred"],
		#theta_ref_pred=ref_data["theta_ref_pred"],
		#ori_phi_ref_pred=ref_data["ori_phi_ref_pred"],
		#hbond_ref_pred=ref_data["hbond_ref_pred"],
		#sidechain_ref_pred=ref_data["sidechain_ref_pred"],
		#phi_ref_pred=ref_data["phi_ref_pred"],
		#psi_ref_pred=ref_data["psi_ref_pred"],
		)
	print("Completed")

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Generate final pred')
	parser.add_argument('--target', type=str, required=True, default="", help='target protein name')
	parser.add_argument('--input_dir', type=str, required=True, default="./input", help='directory containing features')
	parser.add_argument('--out', type=str, default='', help='directory to save prediction')
	parser.add_argument('--dist_ref', type=str, default='./reference/dist_ref', help='directory of dist ref model')
	parser.add_argument('--omega_ref', type=str, default='./reference/omega_ref', help='directory of omega ref model')
	parser.add_argument('--theta_ref', type=str, default='./reference/theta_ref', help='directory of theta ref model')
	parser.add_argument('--ori_phi_ref', type=str, default='./reference/ori_phi_ref', help='directory of ori phi ref model')
	parser.add_argument('--phi_ref', type=str, default='./reference/backbone_phi_ref', help='directory of phi ref model')
	parser.add_argument('--psi_ref', type=str, default='./reference/backbone_psi_ref', help='directory of psi ref model')
	parser.add_argument('--hbond_ref', type=str, default='./reference/hbond_ref', help='directory of hbond ref model')
	parser.add_argument('--sidechain_ref', type=str, default='./reference/sce_ref', help='directory of sidechain ref model')
	args = parser.parse_args()

	target = args.target
	fasta_dir = join(args.input_dir, target)
	out_dir = args.out
	dist_pred_dir = out_dir
	sidechain_pred_dir = fasta_dir
	hbond_pred_dir = fasta_dir
	dist_ref_model_path = args.dist_ref
	omega_ref_model_path = args.omega_ref
	theta_ref_model_path = args.theta_ref
	ori_phi_ref_model_path = args.ori_phi_ref
	phi_ref_model_path = args.phi_ref
	psi_ref_model_path = args.psi_ref
	hbond_ref_model_path = args.hbond_ref
	sidechain_ref_model_path = args.sidechain_ref

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	main(target,fasta_dir,out_dir,dist_pred_dir,sidechain_pred_dir,hbond_pred_dir,dist_ref_model_path, 
		omega_ref_model_path, theta_ref_model_path, ori_phi_ref_model_path, phi_ref_model_path, 
		psi_ref_model_path, hbond_ref_model_path, sidechain_ref_model_path)




		


	
