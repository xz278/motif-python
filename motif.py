"""
motif analysis

"""
import numpy as np
import matplotlib.pyplot as plt
import OpticsClusterArea as OP
from itertools import *
import AutomaticClustering as AutoC
import pickle
from math import radians, cos, sin, asin, sqrt
from sklearn.cluster import DBSCAN


def save_to_file(output_file_name,data):
	with open(output_file_name,'wb') as f:
		pickle.dump(data,f)

def read_from_file(input_file_name):
	with open('input_file_name','rb') as f:
		ret = pickle.load(f)
	return ret

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



class Graph:
	"""
	TO-DO: modify methods for non-interchangable graph

	A graph object stores the information for a store

	Args:
		am: adjacency matrix
		locs: location ids (cluster ids)

	Test example:
	from motif import Graph
	g1 = [[0,0,0],[1,0,1],[1,0,0]]
	g2 = [[0,1,0,1],[1,0,0,0],[0,0,0,0],[1,0,1,0]]
	g3 = [[0,1,1],[1,0,1],[0,0,0]]
	g5 = [[0,1,1],[0,0,1],[0,0,0]]
	locs = [1,2,3]
	g4 = [[0,0,0],[1,0,0],[1,1,0]]
	g = Graph(am = g1)

	"""
	def __init__(self, am = [],locs = []):
		self.am = am
		self.locs = locs
		if len(locs)==0:
			self.isInterch = False
		else:
			self.isInterch = True

		nv = len(am) # number of vertices
		self.outdeg = [0] * nv
		self.indeg = [0] * nv
		for i in range(nv):
			self.outdeg[i] = sum(am[i])
			temp = 0
			for j in range(nv):
				temp = temp + am[j][i]
			self.indeg[i] = temp
		self._reorder()

	def _reorder(self):
		"""
		reorder the adjacency matrix based on descending order of outdegrees and indegrees
		"""
		nv = len(self.am)
		# idx = sorted(range(nv), key = lambda x: self.outdeg[x] * (nv + 1) + self.indeg[x], reverse = True)
		deg = [0] * nv
		for i in range(nv):
			deg[i] = self.outdeg[i] * (nv + 1) + self.indeg[i]
		idx = sorted(range(nv),key = lambda x: deg[x], reverse = True)
		deg.sort(reverse = True)

		# construct permutaion list
		cl = [idx[0]]
		listToCombine = []
		for i in range(1,nv-1):
			if deg[i] != deg[i-1]:
				currPl = self._permuteidx(list(cl)) # current permutaion list
				listToCombine.append(list(currPl))
				cl = [idx[i]]
			else: # if deg[i] == deg[i-1]
				cl.append(idx[i])
		i = nv-1
		if deg[i] == deg[i-1]:
			cl.append(idx[i])
			currPl = self._permuteidx(list(cl))
			listToCombine.append(list(currPl))
		else:
			currPl = self._permuteidx(list(cl))
			listToCombine.append(list(currPl))
			listToCombine.append([[idx[i]]])
		
		permlist = self._get_perm_list(listToCombine)
		nperm = len(permlist)
		self.alt = []
		self.alt_loc = []
		if nperm > 1:
			for i in range(1,nperm):
				temp_alt_am, temp_alt_loc = self._rearrange(permlist[i])
				self.alt.append(temp_alt_am)
				self.alt_loc.append(temp_alt_loc)
		self._change_to(permlist[0])



	def _get_perm_list(self,ls):
		retls = []
		self._get_permlist_helper(ls,0,[],retls)
		return retls

	def _get_permlist_helper(self,ls,pos,currl,retls):
		lslen = len(ls)
		if pos > lslen-1:
			retls.append(list(currl))
			return
		currPmt = ls[pos]
		lc = len(currPmt)
		for i in range(lc):
			self._get_permlist_helper(ls,pos+1,currl + currPmt[i],retls)
				
	def _permuteidx(self,ls):
		retls = []
		self._permuteidx_helper(ls,[],retls)
		return retls


	def _permuteidx_helper(self,remls,currls,retls):
		n = len(remls)
		if n==0:
			retls.append(list(currls))
			return
		for i in range(n):
			temp = remls[i]
			remls.remove(temp)
			currls.append(temp)
			self._permuteidx_helper(remls,currls,retls)
			remls.insert(i,temp)
			currls.remove(temp)


	def _change_to(self,idx):
		"""
		change current adj mat to the one specified by idx
		Args:
			idx: new order list for the adjacency list
		"""
		nv = len(idx)
		temp_am = [row[:] for row in self.am] # alt: use copy() or deepcopy() imported from copy
		temp_outdeg = list(self.outdeg)
		temp_indeg = list(self.indeg)
		temp_loc = list(self.locs)
		for i in range(nv):
			self.outdeg[i] = temp_outdeg[idx[i]]
			self.indeg[i] = temp_indeg[idx[i]]
			if self.isInterch:
				self.locs[i] = temp_loc[idx[i]]
			for j in range(nv):
				self.am[i][j] = temp_am[idx[i]][idx[j]]

	def _rearrange(self,idx):
		"""
		Reorder the adjcency matrix to the specified new adjmat
		Args:
			idx: new order for the adjacency matrix
		Return:
			A new adjacency matrix, without changing current grpah object
		"""
		nv = len(idx)
		newgraph = [row[:] for row in self.am]
		new_loc = [0] * nv
		for i in range(nv):
			new_loc[i] = self.locs[idx[i]]
			for j in range(nv):
				newgraph[i][j] = self.am[idx[i]][idx[j]]
		return newgraph, new_loc

	def ism_to(self,other_graph):
		l1 = len(self.am)
		l2 = len(other_graph.am)
		if l1 != l2:
			return False
		if self.isInterch != other_graph.isInterch:
			return False
		if self.isInterch:
			is_ism = (self.am == other_graph.am) and (self.locs == other_graph.locs)
			if is_ism:
				return True
			l = len(self.alt)
			if l == 0:
				return False
			for i in range(l):
				if (self.alt[i] == other_graph.am) and (self.alt_loc[i] == other_graph.locs):
					return True
			return False
		else:
			is_ism = self.am == other_graph.am
			if is_ism:
				return True
			l = len(self.alt)
			if l == 0:
				return False
			for i in range(l):
				if self.alt[i] == other_graph.am:
					return True
			return False

	def write_motif(self,outputfilename):
		"""
		Write motif info to csv file.
		"""
		nNodes = len(self.am)
		nAlt = len(self.alt)
		if self.isInterch:
			isInterch = 1
		else:
			isInterch = 0
		with open(outputfilename,'w') as f:
			# f.write(str(nNodes) + ',' + str(nAlt) + ',' + str(isInterch) + '\n'); # write mega
			f.write(str(nNodes) + ',' + str(isInterch) + '\n'); # write mega
			# write first adjacency matrix
			for i in range(nNodes):
				templine = ''
				for j in range(nNodes):
					if self.am[i][j] == 1:
						templine += str(j) + ','
				f.write(templine[:-1] + '\n')
			# # write alternative adjacency matrix if exists any
			# if self.isInterch:
			# 	for a in range(nAlt):
			# 		for i in range(nNodes):
			# 			templine = ''
			# 			for j in range(nNodes):
			# 				if self.alt[a][i][j] == 1:
			# 					templine += str(j) + ','
			# 			f.write(templine[:-1] + '\n')
			# write location ids if nodes are interchangable
			if isInterch:
				temploc = ''
				for i in range(nNodes):
					temploc += str(self.locs[i]) + ','
				f.write(temploc[:-1] + '\n')

	@staticmethod
	def read_motif(input_file_name):
		with open(input_file_name,'r') as f:
			lines = f.readlines()
			num_lines = len(lines)
			mega_info = lines[0].split()[0].split(',')
			num_nodes = int(mega_info[0])
			if int(mega_info[1]) == 1:
				isInterch = True
			else:
				isInterch = False
			am = [[0] * num_nodes for _ in range(num_nodes)]
			if isInterch:
				locs = [0] * num_nodes
			else:
				locs = []
			for i in range(num_nodes):
				line_idx = i + 1
				cells = lines[line_idx].split()
				if len(cells)==0:
					continue
				adj_list = cells[0].split(',')
				out_degress = len(adj_list)
				for j in range(out_degress):
					am[i][int(adj_list[j])] = 1
			if isInterch:
				cells = lines[-1].split()[0].split(',')
				for i in range(len(cells)):
					locs[i] = int(cells[i])
			return Graph(am = am, locs = locs)



			


# def haversine(lon1, lat1, lon2, lat2):
def haversine(lat1, lon1, lat2, lon2):
	"""
	calcuate distance between two gps coordinates to meters
	"""
	# convert decimal degrees to radians 
	lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
	# haversine formula 
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
	c = 2 * asin(sqrt(a)) 
	m = 6371000 * c
	return m

class ClusterEngine():
	"""
	This class is used to cluster locations into clusters for motif analysis
	"""
	def __init__(self, data = [], algo = 'dbscan', eps = 10, minpts = 10):
		self.data = data
		self.algo = algo
		self.eps = eps
		self.minpts = minpts
		n = len(data)
		self.labels = [-1] * n
		self.dist_matrix = [[0] * n for _ in range(n)]


	def run(self):
		self._compute_dist_matrix()
		self._run_cluster_algo()	

	def _run_cluster_algo(self):
		"""
		Return clusters, -1 for noise
		"""
		if self.algo == 'dbscan':
			self.labels = DBSCAN(eps=self.eps, min_samples=self.minpts,metric="precomputed").fit_predict(self.dist_matrix)
		if self.algo == 'optics':
			self.labels = self._optics_cluster()


	def _compute_dist_matrix(self):
		"""
		calculate pairwise distance
		"""
		n = len(self.data)
		for i in range(n-1):
			for j in range(i+1,n):
				self.dist_matrix[i][j] = haversine(self.data[i][0],self.data[i][1],self.data[j][0],self.data[j][1])
				self.dist_matrix[j][i] = self.dist_matrix[i][j]

	def write_labels(self,outputfilename):
		with open(outputfilename,'w') as f:
			n = len(self.labels)
			for i in range(n-1):
				f.write(str(self.labels[i]))
				f.write('\n')
			f.write(str(self.labels[-1]))

	def plot_cluster(self):
		X = np.array(self.data)
		n = np.size(X,0)
		fig = plt.figure()
		ax = fig.add_subplot(111)

		ax.plot(X[:,0], X[:,1], 'y.')
		colors = cycle('gmkrcbgrcmk')
		unique_cluster = np.unique(self.labels)
		unique_cluster = unique_cluster[1:]
		clrs = zip(unique_cluster,colors)
		# for i in range(n):
		#     ax.plot(X[i,0],X[i,1], clrs[self.labels[i]][1]+'o', ms=5)
		if len(unique_cluster) != 0:
			for i in range(n):
				c = self.labels[i]
				if c == -1:
					continue
			    ax.plot(X[i,0], X[i,1], clrs[c][1]+'o', ms=5)


		plt.savefig('Graph2.png', dpi=None, facecolor='w', edgecolor='w',
		    orientation='portrait', papertype=None, format=None,
		    transparent=False, bbox_inches=None, pad_inches=0.1)
		plt.show()

	def _optics_cluster(self):
	    x = np.array(self.data)
	    RD, CD, order = OP.optics2(x,self.minpts,self.dist_matrix)
	    RPlot = []
	    RPoints = []
	    num_points = np.size(x,0)
	    for item in order:
	        RPlot.append(RD[item])
	        RPoints.append([x[item][0],x[item][1]])

	    rootNode = AutoC.automaticCluster(RPlot, RPoints)
	    leaves = AutoC.getLeaves(rootNode, [])

	    temp_labels = np.ones((num_points)) * -1
	    cluster_cnt = -1
	    for leaf in leaves:
	        cluster_cnt += 1
	        for v in range(leaf.start, leaf.end):
	            temp_labels[v] = cluster_cnt
	    self.labels = labels
