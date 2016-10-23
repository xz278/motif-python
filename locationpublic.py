import pandas as pd
import numpy as np
import geohash
from collections import Counter
import anvil
from anvil import api
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import folium 
sns.set_context("poster", font_scale=1.4)

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
        
