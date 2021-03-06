"""
motif analysis

"""
import numpy as np
import matplotlib.pyplot as plt
# import OpticsClusterArea as OP
# from itertools import *
# import AutomaticClustering as AutoC
import pickle
from math import radians, cos, sin, asin, sqrt
from sklearn.cluster import DBSCAN
import time
import datetime as dt
# import hdbscan


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
		self.labels = [-1] * n
		self.dist_matrix = [[0] * n for _ in range(n)]


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
		# if show_time:
			print 'Clustering: ' + str(time.time() - start_time) + ' seconds.'


	def _compute_dist_matrix(self, show_time = False):
		"""
		calculate pairwise distance
		"""
		start_time = time.time()
		# n = len(self.data)
		# for i in range(n-1):
		# 	for j in range(i+1,n):
		# 		self.dist_matrix[i][j] = haversine(self.data[i][0],self.data[i][1],self.data[j][0],self.data[j][1])
		# 		self.dist_matrix[j][i] = self.dist_matrix[i][j]
		self.dist_matrix = vectorized_haversine2(self.data)
		if show_time:
			print 'Computing distance matrix: ' + str(time.time() - start_time) + ' seconds.'

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
				else:
					ax.plot(X[i,0], X[i,1], clrs[c][1]+'o', ms=5)
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
			# print cells
			if len(cells) == 1:
				continue
			user_id = cells[0]
			# print user_id
			# print user_id in valid_user
			if has_valid_user and (user_id[:-1] not in valid_user):
				continue

			user_time = dt.datetime.strptime(cells[1][:-5],'%Y-%m-%dT%H:%M:%S')
			# print user_time
			# print cells[2], cells[3]
			user_location = [float(cells[2]),float(cells[3])]
			user_speed = float(cells[6])
			if user_speed > 1: # speed threshold
				continue
			if user_id not in users:
				users[user_id] = User(user_id)
			users[user_id].add(user_time,user_location,user_speed)
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
	def __init__(self, graph, freq = 0):
		self.graph = graph # a Graph object
		self.freq = freq

class User:
	"""
	A User object stores the location data for one user id.
	Could use a tree structure: maybe implement this later
	"""
	def __init__(self, uid):
		self.list_of_MonthlyData = {} # a dictionary later used to store MonthData object
		self.uid = uid
		self.motif = []
		self.data_for_cluster = {}
		self.num_entry = 0

	def add(self,user_time,user_location,user_speed):
		"""
		args:
			user_time: datetime.datetime()
			user_location: list[lat,lon]
			user_speed: float in meters
		"""
		user_month = user_time.month
		if user_month not in self.list_of_MonthlyData:
			self.list_of_MonthlyData[user_month] = MonthlyData(month_value = user_month)
		self.list_of_MonthlyData[user_month].add(user_time,user_location,user_speed)


	def prepare(self):
		M = {'spring':[3,4,5],'summer':[6,7,8],'fall':[9,10,11],'winter':[12,1,2]}
		self.data_for_cluster = {}
		self.data_for_cluster_idx = {}
		self.day_list = {}
		for season in M:
			months_in_curr_season = M[season]
			temp_loc = np.zeros(shape = [0,2], dtype = 'float')
			cluster_start_idx = [] # inclusive
			cluster_end_idx = [] # exclusive
			cluster_day_list = []
			p = 0
			for month in months_in_curr_season:
				if month in self.list_of_MonthlyData:
					month_data = self.list_of_MonthlyData[month]
					add_loc = month_data.get_location()
					if month_data.num_days>0:
						temp_loc = np.append(temp_loc,add_loc,axis = 0)
						cluster_day_list += month_data.list_of_DailyData.values()
						cluster_start_idx.append(p)
						p += month_data.num_days
						cluster_end_idx.append(p)
			if (len(cluster_day_list) >= 10) and (len(temp_loc) > 0):
				self.data_for_cluster[season] = temp_loc
				self.data_for_cluster_idx[season] = [cluster_start_idx,cluster_end_idx]
				self.day_list[season] = cluster_day_list
		if len(self.data_for_cluster) > 0:
			data_size = 0
			for item in self.data_for_cluster.values():
				data_size += len(item)
			self.total_size = data_size
		else:
			self.total_size = 0

	def _load_cluster_engnie(self):
		self.ces = []
		n = len(self.data_for_cluster)
		for i in range(n):
			location_data = self.data_for_cluster.values()[i]
			ce = ClusterEngine(location_data)

	# def total_size(self):
	# 	return self.total_size

	def run_cluster(self):
		n = len(self.data_for_cluster)
		for i in range(n):
			ce = self.ces[i]
			ce.run(show_time = True)
			curr_day_list = self.day_list.values()[i] # day_list for current season
			curr_idx = self.data_for_cluster_idx.values[i]
			l = len(curr_day_list)
			for j in range(l):
				curr_daily_data = curr_day_list[j]
				c = 0
				for p in range(curr_idx[j][0],curr_idx[j][1]):
					curr_daily_data.location_ids[c] = ce.labels[p]
					c += 1

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
			user_size.append(u.total_size)
		sorted_idx = sorted(range(len(users)),key = lambda x:user_size[x],reverse = True)
		new_user_id = [''] * n
		for i in range(n):
			new_user_id[i] = user_id[sorted_idx[i]]
		user_size.sort(reverse = True)
		return new_user_id, user_size


class MonthlyData:
	def __init__(self,month_value = 0):
		self.month_value = month_value
		self.list_of_DailyData = {}
		self.num_days = 0

	def add(self,user_time,user_location,user_speed):
		user_date = str(user_time.date())
		if user_date not in self.list_of_DailyData:
			self.list_of_DailyData[user_date] = DailyData(user_time.date())
			self.num_days += 1
		self.list_of_DailyData[user_date].add(user_time,user_location,user_speed)

	def get_location(self):
		location = np.zeros(shape = [0,2], dtype = 'float')
		removed = []
		for daystr in self.list_of_DailyData:
			curr_daily_data = self.list_of_DailyData[daystr]
			curr_daily_data.prepare()
			if curr_daily_data.is_valid:
				location = np.append(location,curr_daily_data.get_location(), axis = 0)
			else:
				removed.append(daystr)
		for item in removed:
			del(self.list_of_DailyData[item])
			self.num_days -= 1
		return location

class DailyData:
	def __init__(self,data_date):
		"""
		Args:
			list_of_date: list of datetime.date()
			location_data: two-element list storing lat. and lon.
		"""
		self.data_date = data_date
		self.is_sorted = False
		self.data = []
		self.valid_data = []
		self.is_valid = False
		self.num_entry = 0

	def add(self,user_time,user_location,user_speed):
		self.data.append([user_time,user_location[0],user_location[1],user_speed])
		# index:              0    ,   lat:    1    ,    lon:   2    ,     3
		self.is_sorted = False
		self.num_entry += 1

	# def run filer  speed threshold start/end at the same dai ...
	def prepare(self):
		if not self.is_sorted:
			self.data.sort(key = lambda x: x[0])
			self.is_sorted = True
		l = len(self.data)
		self.location_ids = [-1] * l
		self.graphs = []
		num_points_thr = 4 * 16 # duartion threshold is 16 hours minimum
		if l >= num_points_thr:
			self.is_valid = True
		else:
			self.is_valid = False

	def get_location(self):
		location = np.zeros(shape = [self.num_entry,2], dtype = 'float')
		for i in range(self.num_entry):
			location[i,0] = self.data[i][1]
			location[i,1] = self.data[i][2]
		return location


