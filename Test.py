"""
from motif import Graph
g1 = [[0,0,0],[1,0,1],[1,0,0]]
g2 = [[0,1,0,1],[1,0,0,0],[0,0,0,0],[1,0,1,0]]
g3 = [[0,1,1],[1,0,1],[0,0,0]]
g5 = [[0,1,1],[0,0,1],[0,0,0]]
locs = [1,2,3]
g4 = [[0,0,0],[1,0,0],[1,1,0]]
g = Graph(am = g1)

from motif2 import *
g = Graph(am = [[0,1,1],[1,0,1],[0,0,0]], locs = [1,2,3])
m = Motif(g,20)

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
import time
from math import radians, cos, sin, asin, sqrt

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

def vectorized_haversine(X):
	"""
	vectorized haversine function.

	Args:
		X is n-by-2 matrix, storing lat & lon each row
		X should be in the form of numpy.array
	Return:
		D is an n-by-n pairwise distance matrix
	"""

	start_time = time.time()
	# convert decimal degrees to radians
	v_radians = np.vectorize(radians,otypes = [np.float])
	v_sin = np.vectorize(sin,otypes = [np.float])
	v_cos = np.vectorize(cos,otypes = [np.float])
	v_asin = np.vectorize(asin,otypes = [np.float])
	v_sqrt = np.vectorize(sqrt,otypes = [np.float])
	n = np.size(X,0)
	X2 = v_radians(X)
	lat1 = np.tile(X2[:,0],[n,1]).T
	lat2 = np.tile(X2[:,0],[n,1])
	dlat = lat1 - lat2
	lon1 = np.tile(X2[:,1],[n,1]).T
	lon2 = np.tile(X2[:,1],[n,1])
	dlon = lon1 - lon2
	a = v_sin(dlat/2)**2 + v_cos(lat1) * v_cos(lat2) * v_sin(dlon/2)**2
	c = 2 * v_asin(v_sqrt(a))
	m = 6371000 * c
	print 'Vectorized computation: ' + str(time.time() - start_time) + ' seconds.'
	return m

def vectorized_haversine2(X):
	start_time = time.time()
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
	print 'Optimized computation: ' + str(time.time() - start_time) + ' seconds.'
	return D

def compute_dist_matrix(data):
	"""
	calculate pairwise distance
	"""
	start_time = time.time()
	n = np.size(data,0)
	dist_matrix = np.zeros(shape = [n,n], dtype = 'float')
	for i in range(n-1):
		for j in range(i+1,n):
			dist_matrix[i][j] = haversine(data[i][0],data[i][1],data[j][0],data[j][1])
			dist_matrix[j][i] = dist_matrix[i][j]
	print 'Scalar computation: ' + str(time.time() - start_time) + ' seconds.'
	return dist_matrix

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



# from motif2 import *
# valid_user = csv_read('validUser.csv')
# users = load_location_data(filename = 'location_data.csv',valid_user = valid_user)

# u = users.values()[20]
# u.prepare()
# u.run_cluster()
# u._data.values()[0]._ce.plot_cluster()

# User.prepare_all(users)
# sortedid, sortedsize = User.sort_all(users)

# n_lines = len(output_matrix)
# with open('eg.csv','w') as f:
# 	for i in range(n_lines):
# 		w = len(output_matrix[i])
# 		temp_line = output_matrix[i][0]
# 		for j in range(1,w):
# 			temp_line += ',' + str(output_matrix[i][j])
# 		temp_line += '\n'
# 		f.write(temp_line)

# from motif2 import *
# x = read_from_file('eg.csv')
# y = read_from_file('egl.csv')
# ce = ClusterEngine(data = x)
# # ce.run(show_time = True)
# ce.labels = y
# ce.plot_cluster()

# save_to_file('eglabel.csv')

# import datetime as dt
# str_time = '2015/01/28/09/15/42'
# time_format = '%Y/%m/%d/%H/%M/%S'
# dt.datetime.strptime(str_time,time_format)

# def testdef():
from motif2 import *
users = load_valid_location_data()
print(users)
u = users.get(60)
u.run_cluster()
u.run_analysis()
listc = [u,u]
spec,nspec = Motif.combine(listc)

# ------------- Oct 1 ----------------------
from motif2 import *
users = load_valid_location_data()
# U = [41,42,43,44]
U = [32,33,34,35,36,37,38,39,40]
epss = [3,5,8,10,15,20,25,30]
minptss = [3,5,10,15,20,25]
chosen_users = []
for i in U:
	u = users.get(i)
	# u.prepare()
	# print u
	u.load_dist_matrix()
	chosen_users.append(u)

final_results = []
for ceps in epss:
	for cminpts in minptss:
		print 'eps: ' + str(ceps) + '  cminpts: ' + str(cminpts)
		for u in chosen_users:
			print u._uid
			u.set_cluster_para(minpts = ceps, eps = cminpts)
			u.run_cluster()
			u.run_analysis(round_trip = True)
		final_results.append(Motif.combine(chosen_users))

cnt = -1
for ceps in epss:
	for cminpts in minptss:
		cnt += 1
		filename = 'n' + str(ceps) + str(cminpts) + '.csv'
		Motif.write_to_file(filename,final_results[cnt][1])

# ------------- Oct 1 ---------------------
u030_rct@eureka

from motifanalysis import *
users = read_data()
u = users.get(50)
u.usertfile('u50.csv')


# ----------------------------------------------

from motifanalysis import *
# u = User.userffile('u33.csv')
u = User.userffile('u50.csv')
u.form_daily_data()
pr5nt1u
t = u.getd(44)




from motif2 import *
users = load_valid_location_data()
u2 = users.get('u080_rct@eureka')
u2.run_cluster()



t.plot_traj()


from motifanalysis import *
users = read_data()
u = users.get(36)
u.form_daily_data()
u.activity_level()

# send Saeed u33: u030_rct@eureka

from motifanalysis import *
users = read_data('u050_rct@eureka.csv')
u = users.get(0)
u.run_analysis()

d = u.getmd()
d.run_dbscan(minpts = 3)
d.write_labels()

from motifanalysis import *
users = read_data('u004_rct@eureka.csv')
u = users.get(0)
t = u.getmd()
t.run_dbscan()
t.clusters[0].idx_list
t = u.getd(10)
t.run_dbscan()
t.clusters[0].idx_list

from motifanalysis import *
users = read_data('u010_rct@eureka.csv')
u = users.get(0)
u.run_analysis(daily_minpts = 3, global_minpts = 5)
