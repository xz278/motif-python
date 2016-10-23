import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.patches as mpp
import numpy as np
from motif2 import *
from pytz import timezone
from sklearn.cluster import DBSCAN
import scipy.spatial as ss
import pandas as pd

# website used to check actual location:
# http://www.darrinward.com/lat-long/

def timestr(datetime_obj):
	return dt.datetime.strftime(datetime_obj,'%Y/%m/%d/%H/%M/%S')

def read_data(filename = 'valid_location_data.csv'):
	users = Users()
	with open(filename,'r') as f:
		lines = f.readlines()

	for line in lines:
		line = line[:-2]
		cells = line.split(',')
		user_id = cells[0]
		# print cells
		utc = timezone('UTC')
		ny = timezone('America/New_York')
		user_time = dt.datetime.strptime(cells[1],'%Y/%m/%d/%H/%M/%S')
		user_time = utc.localize(user_time)
		user_time = user_time.astimezone(ny)
		user_lat = float(cells[2])
		user_lon = float(cells[3])
		user_speed = float(cells[4])
		users.add(user_id,user_time,user_lat,user_lon,user_speed)
	users.load()
	users.sortuser()
	print(users)
	return users


def cluster_dist_matrix(clusters):
	n = len(clusters)
	D = np.zeros(shape = [n,n], dtype = 'float')
	for i in range(n-1):
		for j in range(i,n):
			D[i,j] = ss.distance.euclidean(clusters[i],clusters[j])
			D[j,i] = D[i,j]
	pd.DataFrame(D)
	return D

class Users:
	def __init__(self):
		self.users = []

	def add(self,uid,time,lat,lon,speed):
		ex = False
		for u in self.users:
			if u.uid == uid:
				u.add(time,lat,lon,speed)
				ex = True
				break
		if not ex:
			self.users.append(User(uid = uid))
			self.users[-1].add(time,lat,lon,speed)

	def load(self):
		for u in self.users:
			u.sortself()
			u.form_daily_data()

	def sortuser(self):
		self.users.sort(key = lambda x: x.get_size(), reverse = True)

	def get(self,i):
		return self.users[i]

	def __str__(self):
		temp = ''
		c = -1
		for u in self.users:
			c += 1
			temp += str(c) + ' - ' + u.uid + ': ' + str(u.get_size()) + '\n'
		return temp




class User:	
	def __init__(self,uid):
		self.uid = uid
		self.data = []
		self.analyzed = False

	
	def write_csv(self):
		filename = self.uid + '.csv'
		with open(filename,'w') as f:
			f.write('time,latitude,longitude\n')
			for i in range(len(self.data)):
				f.write(str(self.data[i][0]) +','+ str(self.data[i][1]) +','+ str(self.data[i][2]) +'\n')




	def add(self,time,lat,lon,speed):
		self.data.append([time,lat,lon,speed])

	# def run_cluster_all(self):
	# 	n = len(self.data)
	# 	x = []
	# 	for i in range(n):
	# 		x.append((self.data[i][1],self.data[i][2]))
	# 	self.ce = ClusterEngine(data = x)
	# 	self.ce.run(show_time = True)
	# 	clusters = np.unique(self.ce.labels)
	# 	sign_places = {}
	# 	for i in range(n):
	# 		c = self.ce.labels[i]
	# 		if c != -1:
	# 			if c not in sign_places:
	# 				sign_places[c] = []
	# 			sign_places[c].append(x[i])
	# 	for key in sign_places:
	# 		temp = sign_places[key]
	# 		sign_places[key] = [np.mean(temp,0),temp]
	# 	self.sign_places = sign_places
	# 	print 'cluster centroids:'
	# 	for key in sign_places:
	# 		ct = sign_places[key][0]
	# 		print "{:10.6f}".format(ct[0]), "{:10.6f}".format(ct[1])


	# def plot_cluster_all(self):
	# 	self.ce.plot_cluster()

	# def export_centroid_all(self):
	# 	filename = uid + '_centroid_all.csv'
	# 	with open(filename,'w') as f:
	# 		for key in self.sign_places:
	# 			ct = self.sign_places[key]
	# 			f.write(str(ct[0]) + ',' + str(ct[1]))
	# 	print 'centroids exported to ' + filename + '.'

	def run_analysis(self,daily_minpts = 3,daily_eps = 10, daily_algo = 'dbscan',global_minpts = 1,global_algo = 'dbscan',global_eps = 10):
		nd = len(self.daily_data)
		self.is_valid = [False] * nd
		th = 24 / 0.25 * 0.8
		temp_data = [] # [centroid]
		temp_info = [] # [daily_date_idx,local_cluster]
		self.nodes_distr = []
		for i in range(len(self.daily_data)):
			dd = self.daily_data[i]
			if len(dd.data) >= th:
				self.is_valid[i] = True
				dd.run_cluster(minpts = daily_minpts, eps = daily_eps, algo = daily_algo)
				# self.nodes_distr.append(dd.num_nodes)
				for k in dd.clusters:
					c = dd.clusters[k]
					temp_data.append(c.centroid)
					temp_info.append([i,k])
			# else:
				# self.nodes_distr.append(0)
		self.global_ce = ClusterEngine(data = temp_data, minpts = global_minpts, eps = global_eps, algo = global_algo)
		self.global_ce.run()
		global_clusters = {}
		for i in range(len(temp_info)):
			global_label = self.global_ce.labels[i]
			d,c = temp_info[i][0],temp_info[i][1]
			dd = self.daily_data[d]
			local_idx = dd.clusters[c].idx_list
			for li in local_idx:
				# print str(dd.start_date), c,len(dd.global_labels), li
				dd.global_labels[li] = global_label
			if global_label != -1:
				if global_label not in global_clusters:
					global_clusters[global_label] = Cluster()
				global_clusters[global_label].add(new_data = temp_data[i])
		for c in global_clusters.values():
			c.gen_centroid()
			c.dis_centroid()
		self.global_clusters = global_clusters
		self.analyzed = True
		print('number of nodes: ' + len(self.global_clusters))
		for i in range(nd):
			if self.is_valid[i]:
				self.nodes_distr.append(self.daily_data[i].get_num_gobal_nodes())
			else:
				self.nodes_distr.append(0)
		maxd = max(self.nodes_distr)
		distr = [0] * (maxd + 1)
		for v in self.nodes_distr:
			distr[v] += 1
		self.nodes_distr = distr
		print('distribution of number of nodes:')
		print('number of nodes      frequency')
		for i in range(maxd):
			tempstr = str(i).rjust(len('number '))
			tempstr += str(self.nodes_distr[i]).rjust(len('of nodes      frequ'))
			print(tempstr)

	def export_global_cluster_centroid(self):
		filename = self.uid + '_cluster.csv'
		if not self.analyzed:
			print('Run analysis() first.')
			return
		with open(filename,'w') as f:
			for gc in self.global_clusters.values():
				f.write(str(gc.centroid[0]) + ',' +str(gc.centroid[1]) + '\n' )
		print('clusters exported to ' + filename + '.')

	def egcc(self):
		self.export_global_cluster_centroid()

	def get_size(self):
		return len(self.data)

	def sortself(self):
		self.data.sort(key = lambda x: x[0])

	def __str__(self):
		temp = self.uid + ',' + str(self.get_size()) + '\n'
		for i in range(self.get_size()):
			temp += timestr(self.data[i][0]) + ',' + str(self.data[i][1]) + ',' + str(self.data[i][2]) + ',' + str(self.data[i][3]) + '\n'
		return temp

	def usertfile(self,filename):
		with open(filename,'w') as f:
			f.write(str(self))

	def time_density(self):
		hours = [0] * 24
		n = len(self.data)
		for i in range(n):
			hours[self.data[i][0].hour] += 1
		# return hours
		outputstr = ' time   frequency\n'
		for i in range(24):
			outputstr += str(i).rjust(len(' tim')) + str(hours[i]).rjust(len('e   freque')) + '\n'
		print(outputstr)


	# def plot_mobility(self):
	# 	n = len(self.data)
	# 	x = []
	# 	y = []
	# 	for i in range(n):
	# 		x.append(self.data[i][1])
	# 		y.append(self.data[i][2])
	# 	plt.scatter(x,y)
	# 	plt.show()

	def plot(self):
		lats = []
		lons = []
		for i in range(len(self.data)):
			lats.append(self.data[i][1])
			lons.append(self.data[i][2])
		plt.scatter(lons,lats)
		plt.show()

	def plot_cluster(self):		
		n = len(self.data)
		d = []
		for i in range(n):
			d.append([self.data[i][1],self.data[i][2]])
		ce = ClusterEngine(data = d)
		ce.run(show_time = True)
		ce.plot_cluster()

	def activity_level(self):
		# haversine(lat1, lon1, lat2, lon2)
		n = len(self.data)
		hours = [0] * 24
		ws = [0] * 24 
		act = [0] * 24
		for i in range(1,n):
			d = haversine(self.data[i-1][1],self.data[i-1][2],self.data[i][1],self.data[i][2])
			h1 = self.data[i-1][0].hour
			h2 = self.data[i][0].hour
			m1 = self.data[i-1][0].minute
			m2 = self.data[i][0].minute
			if h1 == h2:
				hours[h1] += d
				ws[h1] += float(m2 - m1) / 60
			else:
				w1 = float(60 - m1) / 60
				w2 = float(m2) / 60
				hours[h1] += d * w1 / (w1+w2)
				ws[h1] += w1
				hours[h2] += d * w2 / (w1+w2)
				ws[h2] += w2
		for i in range(24):
			act[i] = hours[i] / ws[i]
		outputstr = ''
		for i in range(len(act)):
			outputstr += str(i).rjust(4) + ' : ' + "{:6.1f}".format(act[i]) + ' meters\n'
		print(outputstr)
		# return act

	def getd(self,i):
		return self.daily_data[i]

	def getmd(self):
		ml = 0
		i = -1
		cnt = -1
		for u in self.daily_data:
			cnt += 1
			cl = len(u.data)
			if cl>ml:
				ml = cl
				i = cnt
		return self.daily_data[i]

	def __str__(self):
		temp = 'uid: ' + self.uid + '	\n' 
		n = len(self.daily_data)
		ln = len(str(n))
		temp += 'index'.rjust(ln) + '  number of logs    max/min dist (meter)' + '  total dist (meter)  \n' 
		for i in range(n):
			temp += str(i).rjust(ln) + str(len(self.daily_data[i].data)).rjust(len('  number of')) + ("{:.1f}".format(self.daily_data[i].max_d) + '/' + "{:.1f}".format(self.daily_data[i].min_d)).rjust(len(' logs    max/min')) + "{:.1f}".format(self.daily_data[i].total_dist).rjust(len(' dist (meter)' + '  total d')) + '\n'
		return temp

	def form_daily_data(self):
		n = len(self.data)
		# self.daily_data = []
		temp = {}
		for i in range(n):
			utime = self.data[i][0]
			if utime.hour <= 3:
				tempdatetime = utime - dt.timedelta(days = 1)
			else:
				tempdatetime = utime
			date_str = str(tempdatetime.date())
			if date_str not in temp:
				temp[date_str] = DailyData(tempdatetime.date())
			temp[date_str].add(utime,self.data[i][1],self.data[i][2])
		self.daily_data = sorted(temp.values(),key = lambda x:x.start_date)

	def write_gps_to_csv(self,filename):
		with open(filename,'w') as f:
			for i in range(len(self.data)):
				templine = str(self.data[i][1]) + ',' + str(self.data[i][2]) + '\n'
				f.write(templine)

	@staticmethod
	def userffile(filename):
		with open(filename,'r') as f:
			lines = f.readlines()
		num_lines = len(lines)
		line = lines[0]
		cells = line[:-2].split(',')
		uid = cells[0]
		user = User(uid = uid)
		s = int(cells[1])
		for i in range(1,num_lines):
			line = lines[i]
			cells = line[:-2].split(',')
			time = dt.datetime.strptime(cells[0],'%Y/%m/%d/%H/%M/%S')
			lat = float(cells[1])
			lon = float(cells[2])
			speed = float(cells[3])
			user.add(time,lat,lon,speed)
		user.sortself()
		return user

class DailyData:
	def __init__(self,start_date):
		self.start_date = start_date
		self.time = []
		self.data = np.zeros(shape = [0,2], dtype = 'float')
		self.max_d = 0
		self.min_d = sys.maxint
		self.total_dist = 0
		self.ce_added = False

	def __str__(self):
		temp = ''
		n = len(self.data)
		for i in range(n):
			gps = '(' + str(self.data[i,0]) + ',' + str(self.data[i,1]) + ')'
			temp += str(self.time[i]) + ' : ' + "{:>25}".format(gps)
			if i>0:
				diff_time = (self.time[i] - self.time[i-1]).total_seconds()
				ms = divmod(diff_time,60)
				diff_time_str = ' : ' + "{:6.0f}".format(ms[0]) + ' min ' + "{:2.0f}".format(ms[1]) + ' sec'
				d = haversine(self.data[i-1][0],self.data[i-1][1],self.data[i][0],self.data[i][1])
				s = d / (diff_time/60)
				temp += diff_time_str + ' : ' + "{:8.1f}".format(d) + ' meters : ' + "{:8.1f}".format(s) + ' meter/min'
			temp += '\n'
		return temp

	def run_cluster(self,minpts = -1, eps = -1, algo = ''):
		if not self.ce_added:
			if minpts == -1:
				minpts = 5
			if eps == -1:
				eps = 10
			if len(algo) == 0:
				algo = 'dbscan'
			self.ce = ClusterEngine(algo = algo, data = self.data,eps = eps, minpts = minpts)
			self.ce_added = True
		else:
			if minpts != -1:
				self.ce.set_minpts(minpts)
			if eps != -1:
				self.ce.set_eps(eps)
			if len(algo) != 0:
				self.ce.set_algo(algo)
		self.ce.run()
		self.labels = self.ce.labels
		self.global_labels = [-1] * len(self.data)
		# self.ce.plot_cluster()
		self.clusters = {}

		for i in range(len(self.labels)):
			label = self.labels[i]
			if label != -1:
				if label not in self.clusters:
					self.clusters[label] = Cluster()
				curr_loc = (self.data[i][0],self.data[i][1])
				self.clusters[label].add(new_data = curr_loc,new_idx = i)
		for c in self.clusters.values():
			c.gen_centroid()	
			# c.dis_centroid()
		self.num_nodes = len(self.clusters)

	def get_num_gobal_nodes(self):
		unique_nodes = np.unique(self.global_labels)
		if -1 in unique_nodes:
			l = len(unique_nodes) - 1
		else:
			l = len(unique_nodes)
		return l

	def export_cluster(self):
		filename = str(self.start_date) + '_cluster_centroid'+ '.csv'
		if self.ce_added:
			with open(filename,'w') as f:
				for c in self.clusters.values():
					f.write(str(c.centroid[0]) + ',' + str(c.centroid[1]) + '\n')
			print('cluster centroid exported to ' + filename + '.')
		else:
			print('Run clustering first.')


	# def form_network(self):



	def plot_cluster(self):
		if self.ce_added:
			self.ce.plot_cluster()
		else:
			self.run_dbscan()

	def add(self,time,lat,lon):
		self.time.append(time)
		self.data = np.append(self.data,np.array([[lat,lon]]),axis = 0)
		if len(self.data)>1:
			d = haversine(self.data[-1][0],self.data[-1][1],self.data[-2][0],self.data[-2][1])
			self.total_dist += d
			if d>self.max_d:
				self.max_d = d
			if d<self.min_d:
				self.min_d = d

	def scatter(self):
		plt.scatter(self.data[:,1],self.data[:,0])
		ax = plt.gca()
		ax.set_aspect('equal')
		plt.show()

	# def cluster(self):
	# 	ce = ClusterEngine(data = self.data)
	# 	ce.run(show_time = True)
	# 	ce.plot_cluster()

	def plot_traj(self):
		n = len(self.data)
		# plt.scatter(self.data[:,1],self.data[:,0],color = 'k')
		ax = plt.gca()
		cc = float(n)*0.1
		for i in range(1,n):
			dlon = self.data[i][1] - self.data[i-1][1]
			dlat = self.data[i][0] - self.data[i-1][0]
			ccc = (float(i)+cc)/(float(n)+2*cc)
			c = [ccc,ccc,1]
			# print i,n
			# print c
			# arr = plt.Arrow(self.data[i-1][1], self.data[i-1][0], dlon, dlat, width = 0.000001, headwidth = 0.000005, headlength = 0.00001, fc = c, ec = c)
			arr = mpp.FancyArrow(self.data[i-1][1], self.data[i-1][0], dlon, dlat, width = 0.0000005, length_includes_head = True, head_width = 0.000001, fc = c, ec = c)
			# print self.data[i-1][1], self.data[i-1][0], dlon, dlat
			ax.add_patch(arr)
		# ax.set_axis_bgcolor([0.9,0.9,0.9])
		# plt.scatter(self.data[0,1],self.data[0,0],color = 'y')
		# plt.scatter(self.data[-1,1],self.data[-1,0],color = 'r')
		dlon = self.data[1][1] - self.data[0][1]
		dlat = self.data[1][0] - self.data[0][0]
		arr = mpp.FancyArrow(self.data[0][1], self.data[0][0], dlon, dlat, width = 0.0000005, length_includes_head = True, head_width = 0.000001, fc = 'y', ec = 'y')
		ax.add_patch(arr)
		dlon = self.data[-1][1] - self.data[-2][1]
		dlat = self.data[-1][0] - self.data[-2][0]
		arr = mpp.FancyArrow(self.data[-2][1], self.data[-2][0], dlon, dlat, width = 0.0000005, length_includes_head = True, head_width = 0.000001, fc = 'r', ec = 'r')
		ax.add_patch(arr)
		ax.set_axis_bgcolor([0,0,0])
		ax.autoscale_view(True,True,True)
		ax.set_aspect('equal')
		plt.show()

	def write_to_csv(self,filename):
		with open(filename,'w') as f:
			for i in range(len(self.data)):
				templine = str(self.data[i,0]) + ',' + str(self.data[i,1]) + '\n'
				f.write(templine)

	def wtc(self):
		filename = str(self.start_date) + '.csv'
		self.write_to_csv(filename)
		print('data stored in ' + filename)

	def write_labels(self):
		filename = str(self.start_date) + '_clusterlabels'+ '.csv'
		if self.ce_added:
			with open(filename,'w') as f:
				cl = {}
				for i in range(len(self.data)):
					lb = self.labels[i]
					if lb != -1:
						if lb not in cl:
							cl[lb] = []
						cl[lb].append((self.data[i][0],self.data[i][1]))
				for key in cl:
					f.write('cluster ' + str(key) + ':\n')
					x = cl[key]
					for i in range(len(x)):
						templine = str(x[i][0]) + ',' + str(x[i][1]) + '\n'
						f.write(templine)
					f.write('\n')
			print('clusters exported to ' + filename + '.')
		else:
			print('Run clustering first.')

class Cluster:
	def __init__(self):
		self.data = []
		self.idx_list = []
	
	def add(self,new_data = [0,0], new_idx = 0):
		self.data.append(new_data)
		self.idx_list.append(new_idx)

	def gen_centroid(self):
		self.centroid = np.mean(self.data,0)
	def dis_centroid(self):
		tempstr = str(self.centroid[0]) + ',' + str(self.centroid[1])
		print(tempstr)
	def __str__(self):
		tempstr = ''
		for i in range(len(self.idx_list)):
			tempstr += str(self.idx_list[i]).rjust(3) + ' : ' + str(self.data[i][0] + ',' + str(self.data[i][1]) + '\n')
		return tempstr

def write_to_csv(filecontent, filename):
	x = np.array(filecontent)
	n = np.size(x,0)
	d = np.size(x,1)
	with open(filename,'w') as f:
		for i in range(n):
			templine = ''
			for j in range(d-1):
				templine += str(x[i,j]) + ','
			templine += str(x[i,-1]) + '\n'
			f.write(templine)