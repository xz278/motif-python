if __name__ == '__main__':
	from locationpublic import *
	import os
	import sys
	import eurekacluster
	import math
	from motifanalysis import *
	import motif2

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

	geohash_labels = [-1] * num_points
	q = l.tolist()
	geohash_list = eureka.geo_hash.tolist()
	for i in range(len(geohash_list)):
		current_geohash = geohash_list[i]
		if current_geohash in q:
			geohash_labels[i] = q.index(current_geohash)

	print('No. of unique places: {0}'.format(len(l)))
	x = [geohash.decode(t) for t in l]
	map_view = folium.Map(x[0])

	folium.LatLngPopup().add_to(map_view)
	for p in x:
	    folium.Marker(p).add_to(map_view)


	newpath = './' + uid + '_results'
	if not os.path.exists(newpath):
		os.makedirs(newpath)
	map_view.save(newpath + '/' + uid + '_geohash_cluster_plot.html')


	# Geneate daily network motifs
	nodes = generate_daily_nodes(eureka, 'geo_hash')

	# Only retain days with enough information
	# (i.e., at least 8 non-null intervals)
	nodes_h = [z for z in nodes if len(z[1].dropna()) > 8]
	cols = 15
	rows = 8
	number_graphs = rows * cols
	fig = plt.figure(figsize=(cols*2, rows*2))

	for index, n in enumerate(nodes_h[:number_graphs]):
	    ax = fig.add_subplot(cols, rows, index + 1)
	    g = set(generate_graph(n[1].node)) # removes duplicates
	    nx.draw_spectral(nx.parse_edgelist(g), ax=ax, node_size=30)
	fig.savefig(newpath + '/' + uid + '-graphs-all.pdf')

	eureka_clusters = pd.DataFrame([{'lon': y[1], 'lat': y[0]} for y in x])
	eureka_clusters.to_csv(newpath + '/' + uid + '-geohash-clusters.csv')

	# directed graph
	cols = 15
	rows = 8
	number_graphs = rows * cols
	fig = plt.figure(figsize=(cols*2, rows*2))

	for index, n in enumerate(nodes_h[:number_graphs]):
	    ax = fig.add_subplot(cols, rows, index + 1)
	    g = set(generate_graph(n[1].node))  # removes duplicates
	    dg = nx.DiGraph()
	    
	    # each node in the edge is seperated by whitespace
	    dg.add_edges_from([t.split() for t in g])
	    nx.draw_spectral(dg, ax=ax, node_size=30)
	fig.savefig(newpath + '/' + uid + '-graphs-directed-all.pdf')

	# check frequency
	graphs = {}
	for index, n in enumerate(nodes_h):
	    found = False
	    g = set(generate_graph(n[1].node)) # removes duplicates
	    dg = nx.DiGraph()
	    # each node in the edge is seperated by whitespace
	    dg.add_edges_from([t.split() for t in g]) 
	    dg = nx.freeze(dg)
	    
	    for graph in graphs:
	        if nx.is_isomorphic(graph, dg):
	            # Found an isomorphic graph
	            found = True
	            graphs[graph].append(dg)
	            break
	    if not found:
	        # No isomorphic graph so far, so
	        # create a new one
	        graphs[dg] = [dg]

	# Sort by motif frequency
	d = []
	for k, v in graphs.items():
	    d.append({'g': k, 'total': len(v)})

	d = sorted(d, key=lambda z: z['total'], reverse=True)

	# Generate motif frequency bar plot
	# (limited to top 20 motifs)

	# t = pd.DataFrame(d[:20])
	# t['pct'] = t.total / t.total.sum() * 100
	# fig = plt.figure(figsize=(14, 8))
	# with sns.axes_style('dark'):
	#     ax = fig.add_subplot(111)
	#     t.plot(y='total', kind='bar', ax=ax)
	#     ax.set_xticks([])
	#     ax.legend([])
	#     ax.set_ylabel('Count')
	# fig.savefig(newpath + '/' + uid + '-frequency.pdf')

	# percentage
	t = pd.DataFrame(d[:20])
	t['pct'] = t.total / t.total.sum() * 100
	fig = plt.figure(figsize=(14, 8))
	with sns.axes_style('dark'):
	    ax = fig.add_subplot(111)
	    t.plot(y='pct', kind='bar', ax=ax)
	    ax.set_xticks([])
	    ax.legend([])
	    ax.set_ylabel('Count (%)')
	fig.savefig(newpath + '/' + uid + '-frequency.pdf')

	# Draw motifs sorted by frequency

	fig = plt.figure(figsize=(10, 2))
	for index, g in enumerate(t.g):
	    ax = fig.add_subplot(1, len(t), index + 1)
	    nx.draw_spectral(g, ax=ax, node_size=30)
	fig.savefig(newpath + '/' + uid + '-motif-sorted-by-freq.pdf')



	# First get the top-20 motifs by frequency count
	d = []
	for k, v in graphs.items():
	    d.append({'g': k, 'total': len(v)})
	d = sorted(d, key=lambda z: z['total'], reverse=True)

	# Now sort these motifs by node count

	f = sorted(d[:20], key=lambda z: len(z['g'].nodes()))

	fig = plt.figure(figsize=(10, 2))
	for index, ff in enumerate(f):
	    g = ff['g']
	    ax = fig.add_subplot(1, len(f), index + 1)
	    nx.draw_spectral(g, ax=ax, node_size=30)
	fig.savefig(newpath + '/' + uid + '-motif-sorted-by-nodes.pdf')


	# generate clusters by dbscan
	outputfilename = newpath + '/' + uid + '-dbscan_cluster.csv'
	threshold = math.ceil(num_points * (0.5 / 100))
	dbscan_clusters, ce_clusters = eurekacluster.generate_dbscan_cluster(uid = uid, eps = eps, minpts = threshold, outputfilename = outputfilename)

	geohash_clusters = eureka_clusters.as_matrix(columns = ['lat','lon'])

	start_point = dbscan_clusters[0]
	# print start_point
	cluster_map = folium.Map(location = [start_point[0],start_point[1]])
	folium.LatLngPopup().add_to(cluster_map)

	# add noise points as black circles
	# for i in range(len(self.data)):
	# 	if self.labels[i] == -1:
	# 		loc = [self.data[i][0],self.data[i][1]]
	# 		# print loc
	# 		folium.CircleMarker(location = loc, color = 'black', fill_color = 'black',radius = 3).add_to(cluster_map)

	for c in range(len(dbscan_clusters)):
		centroid = [dbscan_clusters[c][0],dbscan_clusters[c][1]]
		cluster_size = len(ce_clusters[c][1])
		folium.Marker(location = centroid,popup = str(cluster_size), icon = folium.Icon(color = 'red',icon = 'info-sign')).add_to(cluster_map)
	# cluster_map.save(outputfilename)

	for c in range(len(geohash_clusters)):
		centroid = [geohash_clusters[c][0],geohash_clusters[c][1]]
		folium.Marker(location = centroid,icon = folium.Icon(color = 'blue',icon = 'info-sign')).add_to(cluster_map)
	cluster_map.save(newpath + '/' + uid + '-compare-clusters.html')

	print('dbscan cluster distance')
	pd.DataFrame(vectorized_haversine2(dbscan_clusters)).to_csv(newpath + '/' + uid + '-dbscan-cluster-dist.csv')
	print('geohash cluster distance')
	# print(type(geohash_clusters))
	geohash_clusters = geohash_clusters.tolist()
	n = len(geohash_clusters)
	geo_clusters = np.zeros(shape = [n,2], dtype = 'float')
	for i in range(n):
		for j in range(2):
			geo_clusters[i,j] = float(geohash_clusters[i][j])
	pd.DataFrame(vectorized_haversine2(geohash_clusters.tolist())).to_csv(newpath + '/' + uid + '-geohash-cluster-dist.csv')
	motif2.plot_cluster(x = eureka.ix[:,1:], y = geohash_labels, outputfilename = newpath + '/' + uid + '-geohash-cluster-plot.png')