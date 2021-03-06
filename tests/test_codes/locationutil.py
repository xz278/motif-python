import pandas as pd
import numpy as np
import geopy
import geohash
from collections import Counter
import anvil
from anvil import api
import math
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import permutations, repeat
import matplotlib.patches as patches
import pytz
from geopy.distance import vincenty

import networkx as nx
import folium 
import os

def compute_geo_hash(df, lat_c='lat',
                     lon_c='lon', precision=12):
    """
    Computes geohash values.
    
    Parameters
    ----------
    df : DataFrame
    lat_c : str
        Name of column containing latitude values.
        Default is `lat`.
    lon_c : str
        Name of column containing longitude values.
        Default is `lon`.
    precision : int
        Precisoin of generated geohash function.
        Take a look at geohash length and precision
        here: https://en.wikipedia.org/wiki/Geohash
        
        
    Returns
    -------
    l : iterables
        List of geohash values corresponding to the
        rows of values in the given dataframe
    """
    
    
    # get geohash data
    l = []
    for _, row in df.iterrows():
        g = geohash.encode(latitude=row[lat_c],
                           longitude=row[lon_c],
                           precision=precision)
        l.append(g)
        
    return l



def map_geo_hashed_value(l):
    """
    Returns a geohash -> int mapping.
    
    Geohashed values are sorted in ascending order
    and then numerical IDs are assigned.
    
    Parameters
    ----------
    l : iterables
        List of hased values. Might contian
        duplicates.
    
    Returns
    -------
    d : dict
        A dictionary with hashed values as keys
        and numerical ids as values. The dict length
        might be less than the given list if it
        contains duplicates
    """
    
    l = sorted(l)
    return {k: index for index, k in enumerate(l)}

def trim_geo_hash_precision(hashed_values, precision=9):
    """
    Trims geo hash precision.
    
    Parameters
    ----------
    hashed_values : Series
        Contains geohashed values as strings.
        
    precision : int
        The desired precision. If the current
        precision is smaller, then nothing is done.
    """
    
    return hashed_values.map(lambda z: z[:precision])


def filter_out_rare_points(points, threshold_pct=0.5):
    """
    Filters out rare points.
    
    All points with occurrences <= threshold_pct is
    considerd as rare points. All rare points are
    replaced by pd.NaN
    
    Parameters
    ----------
    points : iterables
        Instances of points
        
    threshold_pct : float
        Threshold in percentage for rare points. Any
        point occuring less than given threshold is
        considerd as a rare point. Default is 0.5%
    
    Returns
    -------
    l : list
        List where rare points are marked as pd.NaN.
    """
    
    c = Counter(points)
    total = sum(c.values())
    l = []
    for p in points:
        v = c[p]
        if v/total * 100 <= threshold_pct:
            l.append(np.nan)
        else:
            l.append(p)
    
    return l


def get_primary_location(locations, aggr_f='count',
                         geo_hash_c='geo_hash'):
    """
    Gets the primary location.
    
    Within any given duration there might be a number
    of locations visited by user. This function
    identifies the primary location from a given list
    of locations.
    
    How to define a primary location? In this case,
    we use duration of time spent by a user within
    a given duration.
    
    Parameters
    ----------
    locations : Series
        Series with geo hased values.
        
    aggr_f : str
        Aggregating function. Default is 'count' which
        will result in counting number of rows. If the
        sampling rate is somewhat consistent then just
        counting the number of times a place has been
        recorded is a good approximation of the time
        spent by a user.
        
    geo_hash_c : str
        Column name that contains geo hash values.
        
    Parameters
    ----------
    location : str
        Returns the primary location.
    """
    
    if aggr_f != 'count':
        raise ValueError('Aggregate function {0} is not supported'.format(aggr_f))

    # sorted by size of each group
    return locations.groupby(locations).size().sort_values(ascending=False).index[0]

def generate_daily_nodes(df, hash_c='geo_hash',
                         geo_hash_preicion=None,
                         shift_day_start=None,
                         rare_pt_pct_th=0.5,
                         valid_day_th=8,
                         start_date=None,
                         end_date=None,
                         home=None,
                         **kwargs):
    """
    Parameters
    ----------
    df : DataFrame
        DataFrame with sorted DateTimeIndex.
        
    hash_c : str
        Columns containing geo hash values.
        
    geo_hash_precision : int
        Desired precision of geo hashed values. See
        `trim_geo_hash_precision` for details. If `None`,
        no trimming will happen. Default is None.
        
    shift_day_start : pd.tslib.Timedelta
        The duration by which the start of the day should
        be shifted. For example, in the original paper,
        the day starts at 3:30AM. Provided value should
        be in pd.tslib.TimeDelta (e.g., for 3:30AM start
        of the day, it should be pd.to_timedelta('3.5H')).
        Default is None.
        
    rare_pt_pct_th : float
        Threshold for rare points in percentage. See
        `filter_out_rare_points` for more details. Default
        is 0.5%. If `None`, no filtering happens.
        
    valid_day_th : int
        Minimum number of intervals for a valid day. If
        a day has < valid_day_th intervals, it will be
        discarded. Default is 8.
        
    start_date : TimeStamp
        Start date to generate nodes. If None, the minimum
        day in the DateTimeIndex will be used.
        
    end_date : TimeStamp
        End date to generate nodes. If None, 1 + maximum day
        in the DateTimeIndex will be used.
        
    kwargs
        Arbitrary keyword based arguments passed to
        `generate_nodes` (e.g., time_interval)
        
        
    Returns
    -------
    l : list
        A list containing (date, nodes) pairs. Where nodes
        are represented as a DataFrame returned by `generate_nodes`.
        
        
    Notes
    -----
        The start date and end date is used pd.date_range to
        generate list of days. So, it is important to have
        same timezone information for start_date and end_date
        as in the given DateTimeIndex.
        
    """
    
    df = df.copy().loc[:, [hash_c]]
    
    l = []
    
    if start_date is None:
        d = df.index.min()
        tz = d.tz  # timezone information
        
        start_date = pd.to_datetime(d.date()).tz_localize(tz)
        
    if end_date is None:
        d = df.index.max()
        tz = d.tz  # timezone information
        
        # maximum date + 1
        end_date = pd.to_datetime(d.date()).tz_localize(tz) + pd.to_timedelta('1D')
        
        
    # shifting start of the day
    if shift_day_start is not None:
        start_date = start_date + shift_day_start
        end_date = end_date + shift_day_start
        
        
    if geo_hash_preicion is not None:
        df[hash_c] = trim_geo_hash_precision(df[hash_c], geo_hash_preicion)
    
    # remove rare points
    if rare_pt_pct_th is not None:
        df[hash_c] = filter_out_rare_points(df[hash_c],
                                            rare_pt_pct_th)
    
    # remove NA values (potentially resulting from removing rare points)
    df = df.dropna(subset=[hash_c])
    
    days = pd.date_range(start=start_date, end=end_date, freq='1D')
    for index, rows in enumerate(anvil.utils.get_df_slices(df, days)):
        d = days[index]
        
        nodes = generate_nodes(rows[hash_c], start_time=d, **kwargs)
        
        if (home is not None) and (type(nodes.ix[0,'node']) is not str):
            nodes.ix[0,'node'] = home
        
        if len(nodes) < valid_day_th:
            l.append((d, np.nan))
        else:
            l.append((d, nodes))
            
    return l

def generate_nodes(locations,
                   start_time,
                   end_time=None,
                   time_interval='30Min',
                   valid_interval_th=1):
    """
    Generates motif information from location data.
    
    This function follows the work of Schneider et al.
    (see http://rsif.royalsocietypublishing.org/content/10/84/20130246/)
    
    Parameters
    ----------
    
    locations : Series
        Series with sorted DateTimeIndex and geo hased values
        over a given day.
        
    start_time : pandas.tslib.Timestamp
        Start time to generate time intervals.
        
    end_time : pandas.tslib.Timestamp
        End time for generating time invervals. If None,
        it is set to 24 hours after start_time. Default
        is None.
        
    time_interval : str
        The interval duration. The default is 30 mins
        (resulting in 48 intervals per day). The argument
        is passed to `pd.date_range` as the frequency string.
        
    valid_interval_th : int
        Minimum number of records for a valid interval.
        If an interval has < valid_interval_th rows, it will
        not be considered. Default is 1.
        
    Returns
    -------
    s : Series
        Series with interval order as keys and primary location as
        values. If the interval is not valid (e.g., not having
        sufficient records), it will contain np.nan as value.
    """
    
    if end_time is None:
        end_time = start_time + pd.to_timedelta('1D')

    intervals = pd.date_range(start=start_time,
                              end=end_time, freq=time_interval)
    s = []
    
    for index, t in enumerate(anvil.utils.get_df_slices(locations,
                                             intervals)):
        debug = intervals[index]
            
        if len(t) < valid_interval_th:
            s.append({'node': np.nan, 'time': intervals[index]})
        else:
            s.append({'node': get_primary_location(t), 'time': intervals[index]})
            
    return pd.DataFrame(s)

def resample_location_data(df, time_delta=None,
                          duration_c='duration',
                          time_c='entry_time'):
    """
    Resamples location data.
    
    Location data from iOS client contains duration
    information. This might result in irregular intervals
    between each subsequent records. As we use frequency
    of records in `get_primary_location`, this can result
    in inconsistent results (e.g., given two instances of L1
    with 10 seconds each and one instance of L2 with 25 minutes,
    it will wrongly identify L1 as primary location). Just
    using duration is not adequate to solve this problem as
    duration might overlap multiple intervals. For example,
    in case of 30 minutes intervals, a location record starting
    at 4:15 with duration of 30 minutes needs to be split into
    two records: 4:15–4:30 and 4:30–4:45. Otherwise, the primary
    location might be wrongly attributed to a different location
    (or no location) for 4:30–5:00 interval.
    
    One way of avoiding these issues is to have somewhat uniform
    sampling. In other words, given a desired sampling rate of 1 min,
    if a row has duration of 30 minutes, then spread this record over
    30 rows where each row represent a single minute. In this case,
    the entry times of each record should be updated accrodingly (i.e.,
    entry time of i-th row should be entry_time + i-th minute) and
    the total duration should add up to 30 minutes.
    
    Parameters
    ----------
    
    df : DataFrame
    time_delta : str
        Argument passed to `pd.to_timedelta`. Default is
        1 minute.
    duration_c : str
        Duration column name. The duration must be in seconds.
    time_c : str
        Time column name. Must be in pd.Timestamp format.
        
    Returns
    -------
    DataFrame
        Returns a new dataframe with updated sampling rate.
    """

    l = []
    if time_delta is None:
        time_delta = pd.to_timedelta('60s')
    else:
        time_delta = pd.to_timedelta(time_delta)
    threshold = time_delta.total_seconds()

    for rows in df.iterrows():
        r = rows[1]
        d = r[duration_c]
        entry_time = r[time_c]
        while d > threshold:
            new_r = r.copy()
            new_r[time_c] = entry_time
            new_r[duration_c] = threshold
            l.append(new_r)

            d = d - threshold
            entry_time = entry_time + time_delta

        new_r = r.copy()
        new_r[time_c] = entry_time
        new_r[duration_c] = d
        l.append(new_r)
        
    return pd.DataFrame(l)


def generate_graph(nodes):
    """
    Generate graphs from a given series of nodes.
    
    Parameters
    ----------
    nodes : Series
        An iterable of nodes. It can contain NaN values
        
        
    Returns
    -------
    list
        A list of strings where the edges are seperated by
        whitespace (e.g., ["a b", "b c"]). It can be parsed
        by networkx.parse_edgelist.
        
        
    ToDos
    -----
        1. Only considering non-consecutive nodes (e.g.,
        what to do if seperated by NaN).
    """
    
    nodes = nodes.dropna()
    l = []
    edge_format = "{0} {1}"
    
    for x1, x2 in zip(nodes.shift(), nodes):
        if not(pd.isnull(x1) or pd.isnull(x2)):
            if x1 != x2:
                l.append(edge_format.format(x1, x2))
                
    return l


def get_geo_center(df, lat_c='latitude', lon_c='longitude'):
    """
    Calculates center of given geo points.
    
    To avoid "wrap-around" issues (near international dateline),
    this method follows [1].
    
    Parameters
    ----------
    df: DataFrame
    
    lat_c: str
        Column name containing latitudes. Default is latitude.
    lon_c: str
        Column name containing longitudes. Default is longitude
        
    Returns
    -------
    d: dict
        Keys 'latitude' and 'longitude' indicate
        the coordinates of the center point.

    [1]: https://carto.com/blog/center-of-points/
    """
    
    lat = df[lat_c]
    lon = df[lon_c]
    
    angle = math.pi / 180
    
    center_lat = math.atan2(np.mean(np.sin(lat * angle)), 
                            np.mean(np.cos(lat * angle))) * 180 / math.pi
    
    center_lon = math.atan2(np.mean(np.sin(lon * angle)), 
                            np.mean(np.cos(lon * angle))) * 180 / math.pi
    
    return {'longitude': center_lon, 'latitude': center_lat}


def get_stay_point(df, lat_c='latitude',
                   lon_c='longitude', dist_th=300,
                   time_th='30m'):
    """
    Calculates stay points.
    
    Stay points are determined following [1]. That is,
    if a user stays within a given radius (e.g., dist_th)
    for a given period of time (e.g., time_th), then the
    center of these points is a stay point.
    
    In this implementation, given a set of points p_k (i <= k <= j),
    where these points are sorted ascendingly in terms of occurence
    (i.e., t_k <= t_(k+1)), we take a greedy approach to
    determine if these points result in a stay point. That is, these
    points are considered for stay points if:
    
        1. distance between first (p_i) and last point (p_j) is less
    than dist_th, and
        2. the time between the first and last point is larger than the
    threshold (e.g., t_j - t_i >= time_th).
    
    
    Note that this is a greedy algorithm as we implicitly consider
    the first point as the center of a cluster (with radius of dist_th)
    which might not be the optimal choice (similar case for time duration
    calculation). However, these issues should get resolved in the
    subsequent phase (i.e., by computing stay regions).
    
    Parameters
    ----------
    df: DataFrame
        DataFrame with sorted (ascending) DateTimeIndex
        
    lat_c: str
        Column name with latitude values
    
    lon_c: str
        Column name with longitude values
        
    dist_th: float
        Distance threshold in meters. Default is 300m.
        
    time_th: str or pd.timedelta
        Time threshold that will be parsed by pd.to_timedelta.
        Default is 30 minutes.
        
        
    Returns
    -------
    stay_points : list
        A list corresponding to the rows of df. Each item
        represents the sray point ids (so points within same
        stay points will have same id). If a point is a travel point
        (i.e., does not belog to any stay points), then it will 
        have np.NaN.
        
    
    [1]: https://dl.acm.org/citation.cfm?id=1463477
    """
    
    
    index = 0
    stay_points_c = 0  # total stay points count
    stay_points = []
    
    max_len = len(df)
    time_th = pd.to_timedelta(time_th)
    
    while index < max_len:
        mem_c = 1  # current stay point members: just index
        
        # time diff between current and first members
        time_diff = pd.to_timedelta('0s')
        
        
        p_f = geopy.point.Point(latitude=df.iloc[index][lat_c],
                                longitude=df.iloc[index][lon_c])

        time_f = df.index[index]

        j = index + 1
        while j < max_len:
            p_s = geopy.point.Point(latitude=df.iloc[j][lat_c],
                                    longitude=df.iloc[j][lon_c])
            
            d = geopy.distance.GreatCircleDistance(p_f, p_s).m
            if d <= dist_th:
                mem_c += 1  # new member
                time_diff = df.index[j] - time_f  # update total time spent
                j += 1
            else:
                # spatial constrain is not met
                break
        
        
        # Check if previous points met the time threshold constraint
        if time_diff >= time_th:
            # All these points share same stay point id
            stay_points.extend([stay_points_c] * mem_c)
            stay_points_c += 1
        else:
            # these are not valid stay points
            stay_points.extend([np.nan] * mem_c)

        # points up to j has been considered
        index = j
          
    return stay_points       

def merge_neighboring_grid(geo_hash):
    """
    Merges neighboring grids using a greedy approach.
    
    This function sorts each grid by the number by
    frequency. A grid G is merged with K, if:
        1. G and K is neighboring grid (note that
    the definition of neighboring grid depends on
    the precision of geohash), and
        2. K has higher frequency than G.
        
    The resulting merged grid retains the geohash
    of grid with higher frequency (e.g., in previous
    example, the merged grid will have the same geohash
    of K).
    
    Note that this is a greedy approach, so there is
    no gurantee that this would result in optimal clustering.
    For example, there might be cases in which a more
    cascaded merging would be better (e.g., I -> G -> K, but
    in our approach I and K might not be merged). However, this
    serves to constraint the maximum size of the merged grid
    (e.g., 9 grids) at a time, which has been shown to be a
    better strategy by Zheng et al. (2010):
    http://portal.acm.org/citation.cfm?doid=1772690.1772795
    
    Parameters
    ----------
    geo_hash : Series
        A series with geo hashed values.
        
    
    Returns
    -------
    Series
        A new series with updated geo hash values
        after merging neighboring grids.
    """
    
    
    c =  Counter(geo_hash.dropna())
    d = {}
    
    # sort by frequency
    for z, _ in c.most_common():
        # this check is necessary as we remove
        # items dynamically
        if z in c:
            d[z] = z
            
            # go through the potential merge options
            for n in geohash.neighbors(z):
                if n in c:
                    d[n] = z
                    del c[n]  #  merged with grid z
                    
    return geo_hash.map(d)

def get_stay_region(df, stay_point_c='stay_point',
                    lat_c='latitude', lon_c='longitude',
                    precision=7):
    
    """
    Calculates stay regions.
    
    There are three steps:
    
    1. First compute the center of stay points.
    2. Compute the geohash of stay point centers
    3. Merge neighboring grids.
    
    Parameters
    ----------
    
    df : DataFrame
    
    stay_point_c : str
        Column name with stay point ids. Default
        is `stay_point_id`.
        
    lat_c : str
        Column name with latitude values. Default
        is `latitude`.
        
    lon_c : str
        Column name for longitude values. Default
        is `longitude`.
        
    precision : int
        Geo hash precision. Default is 7.
        
    Returns
    -------
    Series
        Series with stay regions where a stay region
        is defined by a geohash value.
    """
    
    centers = {}
    for k, v in df.groupby(stay_point_c):
        # get stay point centers
        center = get_geo_center(v, lat_c=lat_c, lon_c=lon_c)
        
        # now convert to geo-hash grid
        h = geohash.encode(latitude=center['latitude'],
                           longitude=center['longitude'],
                           precision=precision)
        centers[k] = h
        
    # associate same stay point centers
    # to each record
    stay_points = df[stay_point_c].map(centers)
    
    # now convert the stay points to stay regions
    return merge_neighboring_grid(stay_points)

def add_subplot_axes(ax, rect, axisbg='w'):
    """
    from http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
    example:
    
        fig = plt.figure(figsize=(10,10))
        axes = []
        subpos = [0.2,0.6,0.3,0.3]
        x = np.linspace(-np.pi,np.pi)
        for i in range(4):
            axes.append(fig.add_subplot(2,2,i))
        for axis in axes:
            axis.set_xlim(-np.pi,np.pi)
            axis.set_ylim(-1,3)
            axis.plot(x,np.sin(x))
            subax1 = add_subplot_axes(axis,subpos)
            subax2 = add_subplot_axes(subax1,subpos)
            subax1.plot(x,np.sin(x))
            subax2.plot(x,np.sin(x))
    
    """
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)

    return subax


def get_home_location(df, night_hours=[3,4,5], target_days=range(5)):
    """
    Get the home locations from the the user
    
    Parameters:
    -----------
    df: DataFrame
        location data with datetime index
    
    night_hours: list of integers
        during which hours to collect locations, default is from 3 to 5 am
    
    target_days: list of integers
        during which days to collect home locations, default is weekdays (0 = Monday, ..., 6 = Sunday)
        
    Returns:
    --------
    homes: list of tuples
        (lat, lon) for home locations
    
    home_trajs: list of home trajectories
        a list of list, each sublist contains the gps location for one night in time sequence
    
    """
    df = df.sort_index()
    homes = []
    home_trajs = []
    df2 = df
    df2['date'] = [x.date() for x in df.index]
    df2['hour'] = [x.hour for x in df.index]
    grouped = df2.groupby('date')
    for index,group in grouped:
        if index.weekday() in target_days:
            ext = group.loc[group['hour'].isin(night_hours)]
            if len(ext)>0:
                ext = ext.ix[:,['latitude','longitude']].as_matrix()
                ext = [tuple(x) for x in ext]
                for t in ext:
                    homes.append(t)
                home_trajs.append(ext)
    
    return homes, home_trajs

def approx_home_location(df):
    """
    Calculate the home location based on the assumption that the most frequently visited location
    during 24:00 to 6:00 is the home location.
    
    Parameters:
    -----------
    df: Dataframe
        user location data
    
    Returns:
    --------
    home: String
        home location in geohash value
    """
    
    df2 = df.copy()
    df2['hour'] = [x.hour for x in df2.index]
    df2 = df2.loc[df2['hour'].isin([0,1,2,3,4,5])]
    home = Counter(df2.stay_region.tolist()).most_common()[0][0]
    
    return home



def get_weekday(df):
    """
    retain weekdays
    """
    weekdays = [0,1,2,3,4]
    df2 = df.copy()
    df2['dayofweek'] = [x.dayofweek for x in df.index]

    return df2.loc[df2['dayofweek'].isin(weekdays),:]

def to_csv(uid, df, nodes, time_th, dist_th, geohash_prc, rare_point_th, path='./'):
    """
    Write motif analysis data to two csv files.
    
    Parameters:
    -----------
    df: DataFrame
        user location data, columns = ['latitude', 'longitude', 'time', 'geo_hash', 
                                        'stay_point', 'stay_region', 'dayofweek']
    nodes: List of tuples
        nodes returned from generate_daily_nodes()
    uid: String
        user id
    time_th: String
        time threshold used in stay point detection
    dist_th: Int
        distance threshold used in stay point detection
    rare_point_th: Float
        rare point filter threshold in generate_daily_nodes() function
    geohash_prc: Float
        geohash precision
    """
    
    if path[-1] != '/':
    	path += '/'

    directory = '{}t{}-d{}-r{}-g{}/'.format(path, time_th, dist_th, rare_point_th, geohash_prc)
    uid = uid.split('.')[0]
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    df2 = df.copy()
    df2['tzinfo'] = str(df2.index[0].tzinfo)
    df2 = df2.tz_convert('UTC')
    df2.to_csv(directory+uid+'-data.csv')
    list_df = []
    for node in nodes:
        n = node[1].copy()
        n['timestamp'] = node[0]
        n['tz'] = node[0].tz
        list_df.append(n)
    nodes_df = pd.concat(list_df)
    nodes_df.to_csv(directory+uid+'-nodes.csv')


def from_csv(uid, time_th, dist_th, geohash_prc, rare_point_th, path='./'):
    """
    Restore motif analysis data from csv files created by to_csv()
    
    Parameters:
    -----------
    uid: String
        user id
    time_th: String
        time threshold used in stay point detection
    dist_th: Int
        distance threshold used in stay point detection
    rare_point_th: Float
        rare point filter threshold in generate_daily_nodes() function
    geohash_prc: Float
        geohash precision
    
    Returns:
    --------
    df: Dataframe
        user data containing columns = ['latitude', 'longitude', 'time', 'geo_hash', 
                                        'stay_point', 'stay_region']
    
    nodes: List of tuples
        nodes returned from generate_daily_nodes()
    """
    
    if path[-1] != '/':
    	path += '/'

    directory = '{}t{}-d{}-r{}-g{}/'.format(path,time_th, dist_th, rare_point_th, geohash_prc)
    uid = uid.split('.')[0]
    if not os.path.exists(directory):
        print('Error: file not exists')
        return None, None
        
    df = pd.DataFrame.from_csv(directory+uid+'-data.csv')
    tz_info = df.tzinfo[0]
    df = df.tz_localize('UTC')
    df = df.tz_convert(tz_info)
    df = df.drop('tzinfo', axis=1)
    if 'dayofweek' in df.columns:
        df = df.drop('dayofweek',axis=1)

    df.columns = ['latitude', 'longitude', 'time', 'geo_hash', 
                  'stay_point', 'stay_region']

    list_nodes = pd.DataFrame.from_csv(directory+uid+'-nodes.csv')
    grouped = list_nodes.groupby('timestamp')
    nodes = []
    for index, n in grouped:
        m_timestamp = pd.Timestamp(n.ix[0,'timestamp'], tz=n.ix[0,'tz'], freq='D')
        n = n.drop('timestamp', axis=1)
        n = n.drop('tz', axis=1)
        node = (m_timestamp, n)
        nodes.append(node)
    
    return df, nodes

def read_data(uid):
    """
    Read and process location data.
    
    Parameters:
    uid: String
        file name that stores location data.
        
    Returns:
    df: DataFrame
        location data, set time as index and convert time to local time ('America/New_York')
    """
    
    df = pd.read_csv(uid, usecols=['time','longitude','latitude','geo_hash'])
    df= df.set_index(pd.to_datetime(df.time))
    df = df.tz_localize('UTC')
    df = df.tz_convert('America/New_York')
    df = df.sort_index()
    
    return df
