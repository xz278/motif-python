"""
from motif import Graph
g1 = [[0,0,0],[1,0,1],[1,0,0]]
g2 = [[0,1,0,1],[1,0,0,0],[0,0,0,0],[1,0,1,0]]
g3 = [[0,1,1],[1,0,1],[0,0,0]]
g5 = [[0,1,1],[0,0,1],[0,0,0]]
locs = [1,2,3]
g4 = [[0,0,0],[1,0,0],[1,1,0]]
g = Graph(am = g1)

from motif import *
x = read_data('testd.csv')
ce = ClusterEngine(data = x, algo = 'optics')
ce.run()
ce.plot_cluster()

from motif import *
x = read_data('testd.csv')
ce2 = ClusterEngine(data = x, algo = 'dbscan')
ce2.run()
ce2.plot_cluster()


"""

import numpy as np
import matplotlib.pyplot as plt

def perm(ls):
	retls = []
	helper(ls,[],retls)
	return retls


def helper(remls,currls,retls):
	n = len(remls)
	if n==0:
		retls.append(list(currls))
		return
	for i in range(n):
		temp = remls[i]
		remls.remove(temp)
		currls.append(temp)
		helper(remls,currls,retls)
		remls.insert(i,temp)
		currls.remove(temp)


def get_perm_list(ls):
	"""
	Test example: 
from Test import *
ls = [[[1,2],[2,1]],[[3]],[[4,5],[5,4]]]
r = get_perm_list(ls)
	"""
	retls = []
	get_permlist_helper(ls,0,[],retls)
	return retls

def get_permlist_helper(ls,pos,currl,retls):
	lslen = len(ls)
	if pos == lslen:
		retls.append(list(currl))
		return
	currPmt = ls[pos]
	lc = len(currPmt)
	for i in range(lc):
		# print currl, currPmt
		get_permlist_helper(ls,pos+1,currl + currPmt[i],retls)

def read_data(input_filename):
	"""
	Input: filename for locations
	Output: 2d matrix storing the locations
	"""
	output_matrix = [];
	with open(input_filename,'r') as f:
		lines = f.readlines()
		for line in lines:
			l = line.split()
			fs = l[0].split(',')
			temp = [float(fs[0]),float(fs[1])]
			output_matrix.append(temp)
	return output_matrix

def read_dist_matrix(input_filename):
	output_matrix = []
	with open(input_filename,'r') as f:
		lines = f.readlines()
		for line in lines:
			l = line.split()
			fs = l[0].split(',')
			temp = []
			for t in fs:
				temp.append(float(t))
			output_matrix.append(temp)
	return output_matrix
