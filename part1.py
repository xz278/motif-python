if __name__ == '__main__':
	from locationpublic import *
	import os
	import sys
	import eurekacluster
	import math
	from motifanalysis import *
	import motif2

	print('part1 starts.')
	uid = sys.argv[1]
	config = sys.argv[2] # paras for dbscan	
	paras = pd.read_csv(config,usecols=['minpts','eps'])
	minpts = paras.at[0,'minpts']
	eps = paras.at[0,'eps']

	eureka = pd.read_csv(uid, usecols=['time', 'longitude', 'latitude'])
	num_points = len(eureka)
	# add timezone info
	# eureka = anvil.api.convert_time_zone(df = eureka, column_name = 'time', should_localize = 'America/New_York')
	eureka = anvil.api.convert_time_zone(eureka,'time')
	# generate geo hash
	eureka['geo_hash'] = compute_geo_hash(eureka, lat_c='latitude',lon_c='longitude', precision=7)
	eureka_hash = pd.Series(filter_out_rare_points(eureka.geo_hash))
	l = eureka_hash.dropna().unique()

	print('No. of unique places: {0}'.format(len(l)))

	eureka.to_csv('temp1.csv')
	pd.Series(l).to_csv('temp2.csv')


	print('part1 ends.')	