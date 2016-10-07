import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.patches as mpp
import numpy as np
from motif2 import *


# website used to check actual location:
# http://www.darrinward.com/lat-long/

def timestr(datetime_obj):
	return dt.datetime.strftime(datetime_obj,'%Y/%m/%d/%H/%M/%S')

def read_data():
	filename = 'valid_location_data.csv'
	users = Users()
	with open(filename,'r') as f:
		lines = f.readlines()

	for line in lines:
		line = line[:-2]
		cells = line.split(',')
		user_id = cells[0]
		# print cells
		user_time = dt.datetime.strptime(cells[1],'%Y/%m/%d/%H/%M/%S')
		user_lat = float(cells[2])
		user_lon = float(cells[3])
		user_speed = float(cells[4])
		users.add(user_id,user_time,user_lat,user_lon,user_speed)
	users.load()
	users.sortuser()
	print(users)
	return users

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
	
	def add(self,time,lat,lon,speed):
		self.data.append([time,lat,lon,speed])

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
		return hours

	def plot_mobility(self):
		n = len(self.data)
		x = []
		y = []
		for i in range(n):
			x.append(self.data[i][1])
			y.append(self.data[i][2])
		plt.scatter(x,y)
		plt.show()

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
				ws[h1] += 1
			else:
				w1 = (60 - m1) / 60
				w2 = (m2) / 60
				hours[h1] += d * w1
				ws[h1] += w1
				hours[h2] + d * w2
				ws[h2] += w2
		for i in range(24):
			act[i] = hours[i] / ws[i]
		outputstr = ''
		for i in range(len(act)):
			outputstr += str(i) + ' : ' + "{:6.1f}".format(act[i]) + ' meters\n'
		print outputstr
		# return act

	def getd(self,i):
		return self.daily_data[i]

	def __str__(self):
		temp = 'uid: ' + self.uid + '\n' 
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

	def cluster(self):
		ce = ClusterEngine(data = self.data)
		ce.run(show_time = True)
		ce.plot_cluster()

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