def compute_home_cluster(df,cluster_hash,cluster_gps,save_as):
	"""
	Compute home clusters from a user's geo_hash locations.
	Locations recorded at 3,4,5 a.m. are considered as home.
	Parameters
	----------
	df: DataFrame
		User data.
	cluster_hash: list
		A list of significant locations (in geohash form) generated by geohashing.
	cluster_gps: list
		A list of cooresponding gps coordinates for the clusters in cluster_hash.
	save_as: str
		A string representing the name of the html file generated.
	Returns
	-------
	home_cluster_coord: list
		A list of gps coordinates for home clusters.
	home_cluster_count: list
		A list of frequency of the coorespding clusters in home_cluster_coord.
	"""

	significant_cluster = cluster_hash.tolist()
	timeseries = df.index
	clusterseries = eureka.geo_hash.tolist()
	home_clusters = {}
	n = len(timeseries)
	night_time = [3,4,5]
	for i in range(n):
		h = timeseries[i].hour
		loc = clusterseries[i]
		if (h in night_time) and (loc in significant_cluster):
			if loc not in home_clusters:
				home_clusters[loc] = 0
			home_clusters[loc] += 1

	home_cluster_coord = []
	home_cluster_count = []

	for loc_hash in home_clusters:
		idx = significant_cluster.index(loc_hash)
		home_cluster_coord.append(cluster_gps[idx])
		home_cluster_count.append(home_clusters[loc_hash])


	start_point = home_cluster_coord[home_cluster_count.index(max(home_cluster_count))]
	cluster_map = folium.Map(location = [start_point[0],start_point[1]])
	folium.LatLngPopup().add_to(cluster_map)
	for i in range(len(home_clusters)):
		centroid = home_cluster_coord[i]
		cluster_size = home_cluster_count[i]
		folium.Marker(location = centroid,popup = str(cluster_size), icon = folium.Icon(color = 'red',icon = 'info-sign')).add_to(cluster_map)
	
	cluster_map.save(save_as)
	return home_cluster_coord, home_cluster_count




eureka = pd.read_csv(uid, usecols=['time', 'longitude', 'latitude'])
num_points = len(eureka)
# add timezone info
eureka = anvil.api.convert_time_zone(eureka,'time')
# generate geo hash
eureka['geo_hash'] = compute_geo_hash(eureka, lat_c='latitude',lon_c='longitude', precision=7)
eureka_hash = pd.Series(filter_out_rare_points(eureka.geo_hash))
l = eureka_hash.dropna().unique()
x = [geohash.decode(t) for t in l]
# compute home clusters
home_cluster_corrd, home_cluster_count = compute_home_cluster(eureka,l,x,'/tmp/ureka-viz/ureka-u010-home-cluster.html')