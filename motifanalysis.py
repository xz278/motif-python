import datetime as dt


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
