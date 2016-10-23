# import matplotlib.pyplot as plt
# import OpticsClusterArea as OP
from itertools import *
import matplotlib.pyplot as plt
# import AutomaticClustering as AutoC
# import pickle
from math import radians, cos, sin, asin, sqrt
from sklearn.cluster import DBSCAN
import time
import datetime as dt
# from operator import itemgetter
import sys
import folium
import random
from scipy import spatial
from geopy.distance import vincenty, great_circle
import numpy as np
import pandas as pd
import anvil
# from pytz import timezone1

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
		n = len(self.data)
		self.clusters = {}
		for i in range(n):
			l = self.labels[i]
			if l != -1:
				if l not in self.clusters:
					self.clusters[l] = []
				self.clusters[l].append(self.data[i])
		for c in self.clusters:
			self.clusters[c] = [np.mean(self.clusters[c],0),self.clusters[c]]
		num_cluster = len(self.clusters)
		print('number of clusters: ' + str(num_cluster))

	def export_cluster_centroid_map(self,outputfilename = None):
		if outputfilename is None:
			outputfilename = str(dt.date.today()) + '_cluster_map.html'
		start_point = self.clusters[0][0]
		# print start_point
		cluster_map = folium.Map(location = [start_point[0],start_point[1]])
		folium.LatLngPopup().add_to(cluster_map)
		for c in self.clusters:
			centroid = [self.clusters[c][0][0],self.clusters[c][0][1]]
			print(centroid)
			folium.Marker(location = centroid,icon = folium.Icon(color = 'red',icon = 'info-sign')).add_to(cluster_map)
		cluster_map.save(outputfilename)

	def export_cluster_map(self,outputfilename = None):
		if outputfilename is None:
			outputfilename = str(dt.date.today()) + '_cluster_map.html'
		start_point = self.clusters[0][0]
		# print start_point
		cluster_map = folium.Map(location = [start_point[0],start_point[1]])
		folium.LatLngPopup().add_to(cluster_map)

		# add noise points as black circles
		for i in range(len(self.data)):
			if self.labels[i] == -1:
				loc = [self.data[i][0],self.data[i][1]]
				# print loc
				folium.CircleMarker(location = loc, color = 'black', fill_color = 'black',radius = 3).add_to(cluster_map)

		cluster_colors = []
		for c in self.clusters:
			current_color = "#%06x" % random.randint(0,0xFFFFFF)
			while current_color in cluster_colors:
				current_color = "#%06x" % random.randint(0,0xFFFFFF)
			cluster_colors.append(current_color)
			cluster_points = self.clusters[c][1]
			num_points = len(cluster_points)

			for i in range(num_points):
				current_point = [cluster_points[i][0],cluster_points[i][1]]
				folium.CircleMarker(location = current_point, color = current_color, fill_color = current_color,radius = 5).add_to(cluster_map)
			
			centroid = [self.clusters[c][0][0],self.clusters[c][0][1]]
			# print centroid
			folium.Marker(location = centroid,icon = folium.Icon(color = current_color,icon = 'info-sign')).add_to(cluster_map)



		cluster_map.save(outputfilename)
		print('Map saved to: ')
		print('    ' + outputfilename)


	def _run_cluster_algo(self, show_time = False):
		"""
		Return clusters, -1 for noise
		"""
		start_time = time.time()
		if self.algo == 'dbscan':
			db = DBSCAN(eps=self.eps, min_samples=self.minpts,metric="precomputed").fit(self.dist_matrix)
			self.labels = db.labels_
		if self.algo == 'optics':
			self.labels = self._optics_cluster()
		# if self.algo == 'hdbscan':
		# 	self.labels = hdbscan.HDBSCAN(min_cluster_size = self.minpts).fit_predict(self.dist_matrix)
		if show_time:
			print('Clustering: ' + str(time.time() - start_time) + ' seconds.')


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
			print('Computing distance matrix: ' + str(time.time() - start_time) + ' seconds.')
		self._dist_ready = True

	# def _compute_dist_matrix(self, show_time = False):
	# 	"""
	# 	calculate pairwise distance
	# 	"""
	# 	if self._dist_ready == True:
	# 		return
	# 	start_time = time.time()
	# 	# n = len(self.data)
	# 	# for i in range(n-1):
	# 	# 	for j in range(i+1,n):
	# 	# 		self.dist_matrix[i][j] = haversine(self.data[i][0],self.data[i][1],self.data[j][0],self.data[j][1])
	# 	# 		self.dist_matrix[j][i] = self.dist_matrix[i][j]
	# 	# self.dist_matrix = [[0] * n for _ in range(n)]
	# 	# self.dist_matrix = vectorized_haversine2(self.data)
	# 	d = spatial.distance.pdist(self.data,
	# 							   lambda x,y: vincenty(x,y).m)
	# 	self.dist_matrix = spatial.distance.squareform(d)
	# 	if show_time:
	# 		print 'Computing distance matrix: ' + str(time.time() - start_time) + ' seconds.'
	# 	self._dist_ready = True

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

	def plot_cluster(self,outputfilename = None):
		temp_X = np.array(self.data)  
		n = np.size(temp_X,0)
		X = np.zeros(shape = [n,2], dtype = 'float')
		X[:,0] = temp_X[:,1]
		X[:,1] = temp_X[:,0]
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(X[:,0], X[:,1], 'y.')
		unique_cluster = np.unique(self.labels)
		# unique_cluster = unique_cluster[1:]
		index_to_del = np.where(unique_cluster == -1)[0]
		unique_cluster = np.delete(unique_cluster, index_to_del)
		clrs = np.random.rand(len(unique_cluster),3)
		if len(unique_cluster) != 0:
			for i in range(n):
				c = self.labels[i]
				if c == -1:
					continue
				else:
					# ax.plot(X[i,0], X[i,1], clrs[c][1]+'o', ms=5)
					# print c
					ax.plot(X[i,0], X[i,1],color = clrs[c,:], marker = 'o', ms=5)
		if outputfilename is None:
			outputfilename = dt.datetime.now().strftime('%Y-%m-%d_%H-%M-%S') + '.png'
		plt.savefig(outputfilename, dpi=None, facecolor='w', edgecolor='w',
		    orientation='portrait', papertype=None, format=None,
		    transparent=False, bbox_inches=None, pad_inches=0.1)
		# plt.show()

	def get_parameters(self):
		return self.eps, self.minpts

	def get_eps(self):
		return self.eps

	def get_minpts(self):
		return self.minpts

	def show_parameters(self):
		print('eps = ' + str(self.eps) + ', minpts = ' + str(self.minpts))

	def set_eps(self,new_eps):
		self.eps = new_eps

	def set_minpts(self, new_minpts):
		self.minpts = new_minpts

	def set_algo(self, new_algo):
		self.alog = new_algo

	def set_params(self, new_eps = -1, new_minpts = -1, new_algo = ''):
		if new_eps != -1:
			self.eps = new_eps
		if new_minpts != -1:
			self.minpts = new_minpts
		if len(new_algo) != 0:
			self.algo = new_algo



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
	#     	cluster_cnt += 1
	#         for v in range(leaf.start, leaf.end):
	# 			cluster_idx = order[v]
	# 			temp_labels[cluster_idx] = cluster_cnt
	#     return temp_labels


def generate_dbscan_cluster(uid = None, eps = 10, minpts = 10, outputfilename = None):
	if uid is None:
		print('Inputs must include uid.')
		return
	x = pd.read_csv(uid,usecols=['time','longitude','latitude'])
	x = anvil.api.convert_time_zone(x,'time')
	X = x.as_matrix(columns = ['latitude','longitude'])
	ce = ClusterEngine(data = X, eps = eps, minpts = minpts)
	ce.run(show_time = True)

	ret = []
	if outputfilename is not None:
		n = len(ce.clusters)
		with open(outputfilename, 'w') as f:
			f.write('lat,lon\n')
			for i in range(n):
				centroid = ce.clusters[i][0]
				ret.append(centroid.tolist())
				f.write(str(centroid[0]) + ',' + str(centroid[1]) + '\n')
	else:
		n = len(ce.clusters)
		for i in range(n):
				centroid = ce.clusters[i][0]
				ret.append(centroid.tolist())
	# ce.plot_cluster(outputfilename = outputfilename + '.png')
	return ret, ce.clusters


if __name__ == '__main__':
	import numpy as np
	import pandas as pd

	data = sys.argv[1]
	config = sys.argv[2]
	output = sys.argv[3]
	x = pd.read_csv(data,usecols=['time','longitude','latitude'])
	X = x.as_matrix(columns = ['latitude','longitude'])
	paras = pd.read_csv(config,usecols=['minpts','eps'])
	minpts = paras.at[0,'minpts']
	eps = paras.at[0,'eps']
	print('minpts = ' + minpts + ' eps =', eps)
	ce = ClusterEngine(data = X, eps = eps, minpts = minpts)
	ce.run(show_time = True)
	ce.export_cluster_map(outputfilename = output)


