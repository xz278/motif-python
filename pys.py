def readdata(input_filename):
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

def readdata2(input_filename):
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

from math import radians, cos, sin, asin, sqrt
# def haversine(lon1, lat1, lon2, lat2):
def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
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

from math import radians, cos, sin, asin, sqrt
from pys import readdata
from pys import readdata2
from sklearn.cluster import DBSCAN
e = 10;
m = 10;
X = readdata('testdata.csv')
D = readdata2('D.csv')
db = DBSCAN(eps=e, min_samples=m).fit_predict(X)
db = DBSCAN(eps=e, min_samples=m,metric="precomputed").fit_predict(D)


import hdbscan

clusterer = hdbscan.HDBSCAN(min_cluster_size=10)
cluster_labels = clusterer.fit_predict(x)

hdbscanr = hdbscan.HDBSCAN(min_cluster_size=10).fit_predict(x)