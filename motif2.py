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
import time
import datetime as dt
# import matplotlib.pyplot as plt
from operator import itemgetter
import sys
# import hdbscan


def save_to_file(output_file_name,data):
	with open(output_file_name,'wb') as f:
		pickle.dump(data,f)

def read_from_file(input_file_name):
	with open(input_file_name,'rb') as f:
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
		# print(am)
		# print(locs)
		self.am = am
		self.locs = locs
		if len(locs)==0:
			self.isInterch = False
			self._is_interch = False
		else:
			self.isInterch = True
			self._is_interch = True

		nv = len(am) # number of vertices
		if nv == 1:
			self.alt = []
			self.alt_loc = []
			return
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
		# print(idx)

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
		# print(permlist)
		nperm = len(permlist)
		self.alt = []
		self.alt_loc = []
		if nperm > 1:
			for i in range(1,nperm):
				temp_alt_am, temp_alt_loc = self._rearrange(permlist[i])
				self.alt.append(temp_alt_am)
				self.alt_loc.append(temp_alt_loc)
		self._change_to(permlist[0])

	def get_num_vertices(self):
		return len(self.am)

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

	def set_interch(self,isInterch):
		self._is_interch = isInterch

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
		if self.isInterch:
			new_loc = [0] * nv
		else:
			new_loc = []
		for i in range(nv):
			if self.isInterch:
				new_loc[i] = self.locs[idx[i]]
			for j in range(nv):
				# print 'i is: ' + str(i) + ', j is: ' + str(j)
				# print(self.am)
				newgraph[i][j] = self.am[idx[i]][idx[j]]
		return newgraph, new_loc

	def ism_to(self,other_graph):
		l1 = len(self.am)
		l2 = len(other_graph.am)
		if l1 != l2:
			return False
		if self._is_interch != other_graph._is_interch:
			return False
		if self._is_interch:
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

	def export(self):
		am = self.am
		if self._is_interch:
			locs = self.locs
		else:
			locs = []
		return Graph(am = am, locs = locs)

	def __str__(self):
		if self.isInterch:
			temp = 'Interchangable vertices\n'
		else:
			temp = 'Uninterchangable vertices\n\n'
		temp += 'Adjacency matrix\n'
		n = len(self.am)
		for i in range(n):
			for j in range(n-1):
				temp += str(self.am[i][j]) + ' '
			temp += str(self.am[i][-1]) + '\n'
		if self.isInterch:
			temp += '\ncluster/location id\n'		
			for i in range(n-1):
				temp += str(self.locs[i]) + ' '
			temp += str(self.locs[-1]) + '\n'
		return temp




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


def vectorized_haversine(X):
	"""
	vectorized haversine function.

	Args:
		X is n-by-2 matrix, storing lat & lon each row
		X should be in the form of numpy.array
	Return:
		D is an n-by-n pairwise distance matrix
	"""

	# convert decimal degrees to radians
	v_radians = np.vectorize(radians,otypes = [np.float])
	v_sin = np.vectorize(sin,otypes = [np.float])
	v_cos = np.vectorize(cos,otypes = [np.float])
	v_asin = np.vectorize(asin,otypes = [np.float])
	v_sqrt = np.vectorize(sqrt,otypes = [np.float])
	n = np.size(X,0)
	X2 = v_radians(X)
	lat1 = np.tile(X[:,0],[n,1]).T
	lat2 = np.tile(X[:,0],[n,1])
	dlat = lat1 - lat2
	lon1 = np.tile(X[:,1],[n,1]).T
	lon2 = np.tile(X[:,1],[n,1])
	dlon = lon1 - lon2
	a = v_sin(dlat/2)**2 + v_cos(lat1) * v_cos(lat2) * v_sin(dlon/2)**2
	c = 2 * v_asin(v_sqrt(a))
	m = 6371000 * c
	return m

def vectorized_haversine2(X):
	"""
	vectorized haversine function.

	Args:
		X is n-by-2 matrix, storing lat & lon each row, in np.array
		X should be in the form of numpy.array
	Return:
		D is an n-by-n pairwise distance matrix
	"""

	# convert decimal degrees to radians
	v_radians = np.vectorize(radians,otypes = [np.float])
	v_sin = np.vectorize(sin,otypes = [np.float])
	v_cos = np.vectorize(cos,otypes = [np.float])
	n = np.size(X,0)
	X2 = v_radians(X)
	lat = X2[:,0]
	cos_lat = v_cos(lat)
	lon = X2[:,1]
	D = np.zeros(shape = [n,n], dtype = 'float')
	for i in range(n-1):
		for j in range(i+1,n):
			temp1 = sin((lat[i]-lat[j])/2)**2
			temp2 = cos_lat[i] * cos_lat[j] * sin((lon[i]-lon[j])/2)**2
			temp3 = 6371000 * 2 * asin(sqrt(temp1+temp2))
			D[i,j] = temp3
			D[j,i] = temp3
	return D

class ClusterEngine():
	"""
	This class is used to cluster locations into clusters for motif analysis
	"""
	def __init__(self, data = [], algo = 'dbscan', eps = 10, minpts = 10):
		self.data = np.array(data)
		self.algo = algo
		self.eps = eps
		self.minpts = minpts
		n = len(data)
		self._dist_ready = False
		# self.labels = [-1] * n


	def run(self,show_time = False):
		self._compute_dist_matrix(show_time)
		self._run_cluster_algo(show_time)	

	def _run_cluster_algo(self, show_time = False):
		"""
		Return clusters, -1 for noise
		"""
		start_time = time.time()
		if self.algo == 'dbscan':
			self.labels = DBSCAN(eps=self.eps, min_samples=self.minpts,metric="precomputed").fit_predict(self.dist_matrix)
		# if self.algo == 'optics':
		# 	self.labels = self._optics_cluster()
		# if self.algo == 'hdbscan':
		# 	self.labels = hdbscan.HDBSCAN(min_cluster_size = self.minpts).fit_predict(self.dist_matrix)
		if show_time:
			print 'Clustering: ' + str(time.time() - start_time) + ' seconds.'


	def _compute_dist_matrix(self, show_time = False):
		"""
		calculate pairwise distance
		"""
		if self._dist_ready == True:
			return
		start_time = time.time()
		# n = len(self.data)
		# for i in range(n-1):
		# 	for j in range(i+1,n):
		# 		self.dist_matrix[i][j] = haversine(self.data[i][0],self.data[i][1],self.data[j][0],self.data[j][1])
		# 		self.dist_matrix[j][i] = self.dist_matrix[i][j]
		# self.dist_matrix = [[0] * n for _ in range(n)]
		self.dist_matrix = vectorized_haversine2(self.data)
		if show_time:
			print 'Computing distance matrix: ' + str(time.time() - start_time) + ' seconds.'
		self._dist_ready = True

	def reset_dist_matrix(self):
		self.dist_matrix = []
		self._dist_ready = False

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
		unique_cluster = np.unique(self.labels)
		unique_cluster = unique_cluster[1:]


		## first way to get colosr:
		# colors = cycle('gmkrcbgrcmk')
		# clrs = zip(unique_cluster,colors)
		## alternate way:
		# clrs = []
		# for i in range(len(unique_cluster)):
		# 	temp_clr = np.random.rand(3,1)
		# 	while temp_clr in clrs: 
		# 		temp_clr = np.random.rand(3,1)
		# 	clrs.append(temp_clr)
		clrs = np.random.rand(len(unique_cluster),3)
		# for i in range(n):
		#     ax.plot(X[i,0],X[i,1], clrs[self.labels[i]][1]+'o', ms=5)
		
		if len(unique_cluster) != 0:
			for i in range(n):
				c = self.labels[i]
				if c == -1:
					continue
				else:
					# ax.plot(X[i,0], X[i,1], clrs[c][1]+'o', ms=5)
					ax.plot(X[i,0], X[i,1],color = clrs[c,:], marker = 'o', ms=5)
		plt.savefig('Graph2.png', dpi=None, facecolor='w', edgecolor='w',
		    orientation='portrait', papertype=None, format=None,
		    transparent=False, bbox_inches=None, pad_inches=0.1)
		plt.show()

	def get_parameters(self):
		return self.eps, self.minpts

	def get_eps(self):
		return self.eps

	def get_minpts(self):
		return self.minpts

	def show_parameters(self):
		print 'eps = ' + str(self.eps) + ', minpts = ' + str(self.minpts)

	def set_eps(self,new_eps):
		self.eps = new_eps

	def set_minpts(self, new_minpts):
		self.minpts = new_minpts

	def set_params(self, new_eps = -1, new_minpts = -1):
		if new_eps != -1:
			self.eps = new_eps
		if new_minpts != -1:
			self.minpts = new_minpts

	def load_data(self, new_data):
		self.data = new_data
		n = len(self.data)
		self.labels = [-1] * n
		# self.dist_matrix = [[0] * n for _ in range(n)]

	# def _optics_cluster(self):
	#     x = np.array(self.data)
	#     RD, CD, order = OP.optics2(x,self.minpts,self.dist_matrix)
	#     RPlot = []
	#     RPoints = []
	#     num_points = np.size(x,0)
	#     for item in order:
	#         RPlot.append(RD[item])
	#         RPoints.append([x[item][0],x[item][1]])

	#     rootNode = AutoC.automaticCluster(RPlot, RPoints)
	#     leaves = AutoC.getLeaves(rootNode, [])

	#     temp_labels = np.ones(shape = (num_points), dtype = 'int') * -1
	#     cluster_cnt = -1
	#     for leaf in leaves:
	#         cluster_cnt += 1
	#         for v in range(leaf.start, leaf.end):
	# 			cluster_idx = order[v]
	# 			temp_labels[cluster_idx] = cluster_cnt
	#     return temp_labels


# location_data format:
# uid,time,latitude,longitude,altitude,bearing,speed,travelstate,provider,network_type,accuracy
#  0 , 1  ,   2    ,    3    ,   4    ,   5   ,  6  ,     7     ,   8    ,     9      ,   10
# time format: 2014-02-14T23:13:48.000Z
def csv_read(filename): # return a matrix of strings
	if filename.split('.')[-1] != 'csv':
		filename = filename + '.csv'
	max_columns = 0
	with open(filename,'r') as f:
		lines = f.readlines()
		ret_matrix = []
		for line in lines:
			cells = line[:-2].split(',')
			ret_matrix.append(cells)
			l = len(cells)
			# print l
			if l > max_columns:
				max_columns = l
	if max_columns == 1:
		ret2 = []
		for i in range(len(ret_matrix)):
			ret2.append(ret_matrix[i][0])
		ret_matrix = ret2
	return ret_matrix

def csv_write(filename,output_matrix):
	if filename.split('.')[-1] != 'csv':
		filename = filename + '.csv'
	n_lines = len(output_matrix)
	with open(filename,'w') as f:
		for i in range(n_lines):
			w = len(output_matrix[i])
			temp_line = output_matrix[i][0]
			for j in range(1,w):
				temp_line += ',' + str(output_matrix[i][j])
			temp_line += '\n'
			f.write(temp_line)


def load_location_data(filename,columns = [], valid_user = []):
	data = []
	cnt = 0
	users = {} # dictionary {user_id:User object, ...}
	if len(valid_user) == 0:
		has_valid_user = False
	else:
		has_valid_user = True
	with open(filename,'r') as f:
		lines = f.readlines()
		l_lines = len(lines)
		for i in range(1,l_lines):
			line = lines[i]
			cells = line[:-2].split(',')
			if len(cells) == 1:
				continue
			user_id = cells[0]
			if has_valid_user and (user_id[:-1] not in valid_user):
				continue
			user_time = dt.datetime.strptime(cells[1][:-5],'%Y-%m-%dT%H:%M:%S')
			user_location = [float(cells[2]),float(cells[3])]
			user_speed = float(cells[6])
			if user_speed > 1: # speed threshold
				continue
			if user_id not in users:
				users[user_id] = User(user_id)
			users[user_id].add(user_time,user_location)
	return users

def load_valid_location_data():
	# users = {} # dictionary {user_id:User object, ...}
	# with open('valid_location_data.csv','r') as f:
	# 	lines = f.readlines()
	# 	l_lines = len(lines)
	# 	for i in range(1,l_lines):
	# 		line = lines[i]
	# 		cells = line[:-2].split(',')
	# 		user_id = cells[0]
	# 		# print cells
	# 		user_time = dt.datetime.strptime(cells[1],'%Y/%m/%d/%H/%M/%S')
	# 		user_location = [float(cells[2]),float(cells[3])]
	# 		user_speed = float(cells[4])
	# 		if user_speed > 1: # speed threshold
	# 			continue
	# 		if user_id not in users:
	# 			users[user_id] = User(user_id)
	# 		users[user_id].add(user_time,user_location)
	# return users
	users = Users() # dictionary {user_id:User object, ...}
	with open('valid_location_data.csv','r') as f:
		lines = f.readlines()
		l_lines = len(lines)
		for i in range(1,l_lines):
			line = lines[i]
			cells = line[:-2].split(',')
			user_id = cells[0]
			# print cells
			user_time = dt.datetime.strptime(cells[1],'%Y/%m/%d/%H/%M/%S')
			user_location = [float(cells[2]),float(cells[3])]
			user_speed = float(cells[4])
			if user_speed > 1: # speed threshold
				continue
			users.add_data(user_id,user_time,user_location)
		users.load()
	return users



def write_user_data(filename,users):
	if filename.split('.')[-1] != 'csv':
		filename = filename + '.csv'
	with open(filename,'w') as f:
		f.write('uid,time,lat,lon,speed\n')
		for user in users.values():
			for month_data in user.list_of_MonthlyData.values():
				for daily_data in month_data.list_of_DailyData.values():
					l = len(daily_data.data)
					timestr = dt.datetime.strftime(daily_data.data[0],'%Y-%m-%dT%H:%M:%S')
					for i in range(l):
						f.write(user.uid + ',' + timestr + ',' + daily_data.data[1] + ',' + daily_data.data[2] + ',' + daily_data.data[-1] + '\n')

class Motif:
	"""
	A motif object includes a network/graph, location(if interchangable), frequency
	"""
	def __init__(self, graph, freq = 1):
		self.graph = graph # a Graph object
		self.freq = freq

	def __str__(self):
		"""
		String format:
		number of vertices, whether is interchangable, frequency
		[[adjacency matrix]]
		[locatino/cluster id if applicable]
		"""
		n = len(self.graph.am)
		temp = str(n) + ','
		if self.graph._is_interch:
			temp += str(1) + ','
		else:
			temp += str(0) + ','
		temp += str(self.freq) + '\n'
		for i in range(n):
			for j in range(n-1):
				temp += str(self.graph.am[i][j]) + ','
			temp += str(self.graph.am[i][-1]) + '\n'
		if self.graph._is_interch:
			for i in range(n-1):
				temp += str(self.graph.locs[i]) +','
			temp += str(self.graph.locs[-1]) + '\n'
		return temp

	def equal(self,other_graph):
		return self.graph.ism_to(other_graph)

	def increment(self,inc = 1):
		self.freq += inc
	
	@staticmethod
	def write_to_file(filename,motifs):
		"""
		agrs:
			filename: file to print motifs
			motifs: list of motifs
		"""
		with open(filename,'w') as f:
			for m in motifs:
				f.write(str(m))

	@staticmethod
	def copy(m):
		return Motif(graph = m.graph.export(), freq = m.freq)

	@staticmethod
	def combine(users):
		"""
		args:
			users: a list of user object
		"""
		spec = []
		nspec = []
		for u in users:
			spec_motif, nspec_motif = u.get_motif()
			for m in spec_motif:
				ex = False
				for i in spec:
					if i.equal(m.graph):
						i.increment(m.freq)
						ex = True
						break
				if not ex:
					spec.append(Motif.copy(m))
			for m in nspec_motif:
				ex = False
				for i in nspec:
					if i.equal(m.graph):
						i.increment(m.freq)
						ex = True
						break
				if not ex:
					nspec.append(Motif.copy(m))
		spec.sort(key = lambda x: x.freq, reverse = True)
		nspec.sort(key = lambda x: x.freq, reverse = True)
		return spec, nspec

class Users:
	def __init__(self):
		self._data = {}
		self._is_sorted = False

	def add_user(self,user_id):
		if user_id not in self._data:
			self._data[user_id] = User(user_id)
	def add_data(self,user_id,user_time,user_location):
		if user_id not in self._data:
			self._data[user_id] = User(user_id)
		self._data[user_id].add(user_time,user_location)

	def load(self):
		if self._is_sorted == True:
			return
		User.prepare_all(self._data)
		self._sorted_id, self._sorted_size = User.sort_all(self._data)
		self._is_sorted = True
	
	def cluster_all(self):
		for user in self._data.values():
			user.run_cluster()

	def run_analysis(self,round_trip = False):
		for user in self._data.values():
			user.run_analysis(round_trip = round_trip)

	def __str__(self):
		ostr = ''
		n = len(self._data)
		if self._is_sorted == False:
			self.load()
		ostr += 'In sorted order:\n'
		for i in range(n):
			ostr += str(i) + ' - ' + self._sorted_id[i] + ': ' + str(self._sorted_size[i]) + '\n'
		return ostr

	def get(self,u):
		if type(u) is str:
			if u in self._data:
				return self._data[u]
			else:
				print 'User id not found.\n'
				return 'User id not found.\n'
		if type(u) is int:
			if self._is_sorted:
				return self._data[self._sorted_id[u]]
			else:
				print 'Data not in sorted order. User users.load() first.\n'
				return 'Data not in sorted order. User users.load() first.\n'
		print 'argument has to be either type int or str. \n'
		return 'argument has to be either type int or str. \n'

class User:
	S2M = {'spring':[3,4,5],'summer':[6,7,8],'fall':[9,10,11],'winter':[12,1,2]}
	M2S = {3:'spring',4:'spring',5:'spring', 6:'summer',7:'summer',8:'summer', 9:'fall',10:'fall',11:'fall', 12:'winter',1:'winter',2:'winter'}
	def __init__(self,uid):
		self._uid = uid
		self._data = {}
	def add(self,user_time,user_location):
		cseason = User.M2S[user_time.month]
		if cseason not in self._data:
			self._data[cseason] = DataEngine(cseason)
		self._data[cseason].add(user_time,user_location)

	def prepare(self):
		for data_engine in self._data.values():
			data_engine.prepare()
	def run_cluster(self):
		for data_engine in self._data.values():
			data_engine.run_cluster()

	def get_size(self):
		total_size = 0
		for data_engine in self._data.values():
			total_size += data_engine.get_size()
		return total_size

	def run_analysis(self, round_trip = False):
		for data_engine in self._data.values():
			data_engine.run_analysis(round_trip = round_trip)
		spec = []
		nspec = []
		for des in self._data.values():
			if des.is_sufficient():
				for m in des._nspec_motif:
					ex = False
					for i in nspec:
						if i.equal(m.graph):
							i.increment(m.freq)
							ex = True
							break
					if not ex:
						nspec.append(Motif.copy(m))
				for m in des._spec_motif:
					ex = False
					for i in spec:
						if i.equal(m.graph):
							i.increment(m.freq)
							ex = True
							break
					if not ex:
						spec.append(Motif.copy(m))
		spec.sort(key = lambda x: x.freq, reverse = True)
		nspec.sort(key = lambda x: x.freq, reverse = True)
		self._spec_motif = spec
		self._nspec_motif = nspec


	def get_motif(self):
		return self._spec_motif, self._nspec_motif

	def __str__(self):
		temp = self._uid + '\n'
		cnt  = -1
		for de in self._data:
			cnt += 1
			temp += str(cnt) + ' - ' + de + ': ' + str(len(self._data[de]._ce.data)) + '\n'
		return temp

	def get(self,u):
		if type(u) is str:
			if u in self._data:
				return self._data[u]
			else:
				print 'User id not found.\n'
				return 'User id not found.\n'
		if type(u) is int:
			if u < len(self._data):
				return self._data.values()[u]
			else:
				return 'Index out of bound.\n'
		return 'argument has to be either type int or str. \n'

	def get_cluster_size_distribution(self):
		max_num = 0
		for de in self._data.values():
			if not de.is_sufficient():
				continue
			t = de.get_daily_cluster_num()
			if max(t) > max_num:
				max_num = t
		distribution = [0] * max_num
		for de in self._data.values():
			if not de.is_sufficient():
				continue
			for daily_data in de._data.values():
				distribution[daily_data.get_num_cluster()] += 1

	@staticmethod
	def prepare_all(users):
		for u in users.values():
			u.prepare()
	@staticmethod
	def sort_all(users):
		user_id = users.keys()
		user_obj = users.values()
		user_size = []
		n = len(user_id)
		for u in user_obj:
			user_size.append(u.get_size())
		sorted_idx = sorted(range(len(users)),key = lambda x:user_size[x],reverse = True)
		new_user_id = [''] * n
		for i in range(n):
			new_user_id[i] = user_id[sorted_idx[i]]
		user_size.sort(reverse = True)
		return new_user_id, user_size

	@staticmethod
	def print_users_size(user_id,user_size):
		n_user = len(user_id)
		for i in range(n_user):
			print str(i+1) + ' - ' + user_id[i] + ': ' + str(user_size[i])

class DataEngine:
	def __init__(self,season):
		self._season = season
		self._data = {}
		self._ce = ClusterEngine()
		self._is_sufficient = False

	def add(self,user_time,user_location):
		date_str = str(user_time.date())
		if date_str not in self._data:
			self._data[date_str] = DailyData(user_time.date())
		self._data[date_str].add(user_time,user_location)

	def get_ce(self):
		return self._ce

	def plot_data(self):
		plt.scatter(self._ce.data[:,0],self._ce.data[:,1])
		plt.show()

	def prepare(self):
		cluster_data = np.zeros(shape = [0,2], dtype = 'float')
		p = 0
		for daily_data in self._data.values():
			daily_data.prepare()
			if daily_data.is_valid():
				cluster_data = np.append(cluster_data,daily_data.get_cluster_data(),axis = 0)
				daily_data.set_start(p)
				p += daily_data.get_size()
				daily_data.set_end(p)
		if len(cluster_data) <= 10:
			self._is_sufficient = False
		else:
			self._is_sufficient = True
			self._ce.load_data(cluster_data)

	def is_sufficient(self):
		return self._is_sufficient

	def run_cluster(self):
		if not self._is_sufficient:
			print 'not enough data\n'
		else:
			self._ce.run(show_time = True)
			for daily_data in self._data.values():
				if daily_data.is_valid():
					daily_data.set_cluster(self._ce.labels[daily_data.get_start():daily_data.get_end()])
	def get_size(self):
		return len(self._ce.data)

	def reset_dist_matrix(self):
		self._ce.reset_dist_matrix()

	def run_analysis(self, round_trip = False):
		self._nspec_motif = []
		self._spec_motif = []
		for daily_data in self._data.values():
			if daily_data.is_valid():
				if round_trip and not daily_data.is_round_trip():
					continue
				daily_data.gen_network()
				daily_data._graph.set_interch(False)
				nspec_exist = False
				for motif in self._nspec_motif:
					if motif.equal(daily_data._graph):
						motif.increment()
						nspec_exist = True
						break
				if not nspec_exist:
					self._nspec_motif.append(Motif(graph = daily_data._graph.export()))
				daily_data._graph.set_interch(True)
				spec_exist = False
				for motif in self._spec_motif:
					if motif.equal(daily_data._graph):
						motif.increment()
						spec_exist = True
						break
				if not spec_exist:
					self._spec_motif.append(Motif(graph = daily_data._graph.export()))
				self._nspec_motif.sort(key = lambda x: x.freq, reverse = True)
				self._spec_motif.sort(key = lambda x: x.freq, reverse = True)

		def get_motif(self):
			return self._nspec_motif, self._spec_motif
		def get_daily_cluster_num(self):
			ret = []
			for daily_data in self._data.values():
				ret.append(daily_data.get_num_cluster)
			return ret
				
	@staticmethod
	def _add_network(motif_ls, freq_ls, new_network):
		n = len(motif_ls)
		motif_exist = False
		for i in range(n):
			if new_network.is_ism(motif_ls[i]):
				freq_ls[i] += 1
				motif_exist = True
				break
		if not motif_exist:
			freq_ls.append(new_network)


class DailyData:
	def __init__(self,new_date):
		self._date = new_date
		self._is_valid = False
		self._data = [] # [datetime, lat, lon]
		self._cluster = []
		self._is_sorted = False
		self._start = -1
		self._end = -1
		self._graph = []

	def get_size(self):
		return len(self._data)

	def add(self,user_time,user_location):
		self._data.append([user_time,user_location[0],user_location[1]])
		self._is_sorted = False

	def prepare(self):
		if not self._is_sorted:
			self._data.sort(key = lambda x: x[0])
			self._is_sorted = True
		l = len(self._data)
		self._cluster = [-1] * l
		num_points_thr = 4 * 16 # duartion threshold is 16 hours minimum
		if l >= num_points_thr:
			self._is_valid = True
		else:
			self._is_valid = False

	def gen_network(self):
		if self.is_valid():
			clusters = np.unique(self._cluster)
			if -1 in clusters:
				clusters = np.delete(clusters,np.where(clusters==-1))
			clusters = clusters.tolist()
			num_vertices = len(clusters)
			tnetwork = [[0] * num_vertices for _ in range(num_vertices)]
			num_entry = self.get_size()
			prev_v = -1
			for i in range(num_entry):
				if self._cluster[i] == -1:
					continue
				if prev_v == -1:
					prev_v = clusters.index(self._cluster[i])
					continue
				curr_v = clusters.index(self._cluster[i])
				if curr_v == prev_v:
					continue
				else:
					tnetwork[prev_v][curr_v] = 1
				prev_v = curr_v
			self._graph = Graph(am = tnetwork, locs = clusters)

	def is_valid(self):
		return self._is_valid

	def get_daily_networks(self):
		return self._graphs

	def get_cluster_data(self):
		temp_data = []
		for i in range(len(self._data)):
			loc = [self._data[i][1],self._data[i][2]]
			temp_data.append(loc)
		return np.array(temp_data,dtype = 'float')

	def get_start(self):
		return self._start

	def get_end(self):
		return self._end

	def set_start(self,start):
		self._start = start

	def set_end(self,end):
		self._end = end

	def set_cluster(self,cluster):
		self._cluster = cluster

	def is_round_trip(self):
		i = 0
		l = len(self._data)
		while i<l and self._cluster[i]==-1:
			i +=1
		j = -1
		while j>=-l and self._cluster[j]==-1:
			j -= 1
		if self._cluster[i]==self._cluster[j]:
			return True
		else:
			return False

	def get_num_cluster(self):
		# clusters = np.unique(self._cluster)
		# if -1 in clusters:
		# 	clusters = np.delete(clusters,np.where(clusters==-1))
		# return len(clusters)
		return self._graph.get_num_vertices()