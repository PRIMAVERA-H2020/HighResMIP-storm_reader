"""
Provides example functions that are useful for assessing model 
tropical storms.


"""
import numpy
import datetime, cftime
import calendar

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import shapely.geometry as sgeom

import iris
import iris.coord_systems as icoord_systems
import iris.coords as icoords


#: Lat/lon locations for each ocean basin for mapping. If set to 
#: None then it returns a global map
MAP_REGION = {'na': (-105, 0, 60, 0),
              'ep': (-170, -80, 40, 0),
              'wp': (-265, -180, 50, 0),
              'ni': (-310, -260, 30, 0),
              'si': (-340, -260, 0, -40),
              'au': (-270, -195, 0, -40),
              'sp': (-200, -100, 0, -40),
              'sa': ( -90, 0, 0, -40),
              'nh': (-360, 30, 70, 0),
              'sh': (-360, 30, 0, -90),
              None: (-360, 0, 90, -90)
              }

#: Lat/lon locations of tracking regions for each ocean basin. If None
#: then returns a region for the whole globe
TRACKING_REGION = {'na': ([-75, -20, -20, -80, -80, -100, -100, -75, -75], [0, 0, 60, 60, 40, 40, 20, 6, 0]),
                   'ep': ([-140, -75, -75, -100, -100, -140, -140], [0, 0, 6, 20, 30, 30, 0]),
                   #'wp': ([-260, -180, -180, -260, -260], [0, 0, 60, 60, 0]),
                   'wp': ([-260, -180, -180, -260, -260], [0, 0, 30, 30, 0]),
                   'cp': ([-180, -140, -140, -180, -180], [0, 0, 50, 50, 0]),
                   'ni': ([-320, -260, -260, -320, -320], [0, 0, 30, 30, 0]),
                   'si': ([-330, -270, -270, -330, -330], [-40, -40, 0, 0, -40]),
                   'au': ([-270, -200, -200, -270, -270], [-40, -40, 0, 0, -40]),
                   'sp': ([-200, -120, -120, -200, -200], [-40, -40, 0, 0, -40]),
                   'sa': ([-90, 0, 0, -90, -90], [-40, -40, 0, 0, -40]),
#                   'nh': ([-360, 0, 0, -360, -360],[0, 0, 90, 90 ,0]),
                   'nh': ([-359.9, 0, 0, -359.9, -359.9],[0, 0, 90, 90 ,0]),
#                   'sh': ([-360, 0, 0, -360, -360],[-90, -90, 0, 0 ,-90]),
                   'sh': ([-359.9, 0, 0, -359.9, -359.9],[-90, -90, 0, 0 ,-90]),
                   'mdr': ([-80, -20, -20, -80, -80], [10, 10, 20, 20, 10]),
                   None: ([-360, 0, 0, -360, -360],[-90, -90, 90, 90 ,-90])
                   }    

#: Corresponding full basin names for each abbreviation
BASIN_NAME = {'na': 'North Atlantic',
              'ep': 'Eastern Pacific',
              'wp': 'Western Pacific',
              'cp': 'Central Pacific',
              'ni': 'North Indian Ocean',
              'si': 'Southwest Indian Ocean',
              'au': 'Australian Region',
              'sp': 'South Pacific',
              'nh': 'Northern Hemisphere',
              'sh': 'Southern Hemisphere',
              None: 'Global'
              }

#: Corresponding month name for a given integer value
NUM_TO_MONTH = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                7: 'Jul', 8: 'Aug', 9: 'Sep',10: 'Oct',11: 'Nov',12: 'Dec'}


def _get_time_range(year, months, calendar = 'proleptic_gregorian'):
    """ 
    Creates a start and end date (a datetime.date timestamp) for a 
    given year and a list of months. If the list of months overlaps into 
    the following year (for example [11,12,1,2,3,4]) then the end date 
    adds 1 to the original year 
    
    """
    #start_date = datetime.datetime(year, months[0], 1)
    #print 'get_time_range ',start_date
    start_date = cftime.datetime(year, months[0], 1)
    cdftime = cftime.utime('hours since 1950-01-01 00:00:00', calendar = calendar)
    t = cdftime.date2num(start_date)
    t1 = cdftime.num2date(t)
    
    end_year = year
    end_month = months[-1]+1
    if months[-1]+1 < months[0] or months[-1]+1 == 13 or len(months) >= 12:
        end_year = year+1
    if months[-1]+1 == 13:
        end_month = 1
    #end_date = datetime.datetime(end_year, end_month, 1)
    end_date = cftime.datetime(end_year, end_month, 1)
    t = cdftime.date2num(end_date)
    t2 = cdftime.num2date(t)
    return t1, t2
                
        
def _storms_in_time_range(storms, year, months):
    """Returns a generator of storms that formed during the desired time period """
    calendar = 'proleptic_gregorian'
    for storm in storms[:1]:
        # derive the calendar from the storm object, and then pass this to ensure that the start/end period has the same calendar for comparison
        cal_type = str(type(storm.genesis_date()))
        cal = cal_type.split('.')[-1][8:]
        if '360' in cal:
            calendar = '360_day'
        elif '365' in cal:
            calendar = '365_day'
        elif 'noleap' in cal or 'NoLeap' in cal:
            calendar = 'noleap'
        else:
            calendar = 'proleptic_gregorian'

    start_date, end_date = _get_time_range(year, months, calendar = calendar)

    for storm in storms:        
        #print start_date, end_date, storm.genesis_date()
        if (storm.genesis_date() >= start_date) and (storm.genesis_date() < end_date):
            yield storm


def load_map(basin=None):
    """ Produces map for desired ocean basins for plotting. """ 
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-160))
    if basin == None:
        ax.set_global()
    else:
        ax.set_extent(MAP_REGION.get(basin))
    resolution ='50m' # use '10m' for fine scale and '110m' for  coarse scale (default)    
    return ax
    
    
def _basin_polygon(basin, project=True):
    """ 
    Returns a polygon of the tracking region for a particular 
    ocean basin. i.e. a storm must pass through this region 
    in order to be retined. For example, if basin is set to 
    'au' then storms for the Australian region must pass
    through the area defined by -270 to -200W, 0 to -40S.
    
    """
    rbox = sgeom.Polygon(list(zip(*TRACKING_REGION.get(basin))))
    if project: 
        rbox = ccrs.PlateCarree().project_geometry(rbox, ccrs.PlateCarree())
    return rbox
    
    
def _storm_in_basin(storm, basin):
    """ Returns True if a storm track intersects a defined ocean basin """
    rbox = _basin_polygon(basin)   
    lons, lats = list(zip(*[(ob.lon, ob.lat) for ob in storm.obs]))
    try:
        track = sgeom.LineString(list(zip(lons, lats)))       
        projected_track = ccrs.PlateCarree().project_geometry(track, ccrs.Geodetic())
        if rbox.intersects(projected_track):
            return True
        return False
    except:
        print('failed storm in basin, ',lons, lats)
        return False

def _storm_genesis_in_basin(storm, basin):
    """ Returns True if the maximum intensity of the storm occurred
    in desired ocean basin. """
    rbox = _basin_polygon(basin)  
    if 'obs_at_genesis' in dir(storm):
        xy = ccrs.PlateCarree().transform_point(storm.obs_at_genesis().lon, storm.obs_at_genesis().lat, ccrs.Geodetic())
    else:
        xy = ccrs.PlateCarree().transform_point(storm.obs_at_genesis().lon, storm.obs_at_genesis().lat, ccrs.Geodetic())
    point = sgeom.Point(xy[0], xy[1])
    if point.within(rbox):
        return True
    return False

def _get_genesis_months(storms, years, basin):
    """ 
    Returns genesis month of all storms that formed within a 
    given set of years 
    
    """
    genesis_months = []
    for storm in storms:
        if (storm.genesis_date().year in years) and _storm_in_basin(storm, basin):
            genesis_months.append(storm.genesis_date().month)
    return genesis_months
            
            
def _month_names(months):
    """ Returns list of month names for a given set of integer values """
    names = []
    for month in months:
        names.append(NUM_TO_MONTH.get(month))
    return names
             
def _get_time_period(years, months):
    """ 
    Returns string of time period for a given set of 
    years and months. E.g. months [6,7,8,9] and years 
    numpy.arange(1989,2003) would return a string 
    'June-September 1989-2002'. Note: years and 
    months must be lists or arrays.
    
    
    """    
    start_mon = calendar.month_name[months[0]]
    end_mon = calendar.month_name[months[::-1][0]]
    start_yr = str(years.min())
    end_yr = str(years.max())
    if start_yr == end_yr:
        return '%s-%s %s' % (start_mon, end_mon, start_yr)
    else:
        return '%s-%s %s-%s' % (start_mon, end_mon, start_yr, end_yr)
    
    
def _cube_data(data):
    """Returns a cube given a list of lat lon information."""
    cube = iris.cube.Cube(data)
    lat_lon_coord_system = icoord_systems.GeogCS(6371229)
    
    step = 4.0
    start = step/2
    count = 90
    pts = start + numpy.arange(count, dtype=numpy.float32) * step
    lon_coord = icoords.DimCoord(pts, standard_name='longitude', units='degrees', 
                                 coord_system = lat_lon_coord_system, circular=True)
    lon_coord.guess_bounds()
    
    start = -90
    step = 4.0
    count = 45
    pts = start + numpy.arange(count, dtype=numpy.float32) * step
    lat_coord = icoords.DimCoord(pts, standard_name='latitude', units='degrees', 
                                 coord_system = lat_lon_coord_system)
    lat_coord.guess_bounds()
    
    cube.add_dim_coord(lat_coord, 0)
    cube.add_dim_coord(lon_coord, 1)
    return cube 


def _binned_cube(lats, lons):
    """ Returns a cube (or 2D histogram) of lat/lons locations. """   
    data = numpy.zeros(shape=(45,90))
    binned_cube = _cube_data(data)
    xs, ys = binned_cube.coord('longitude').contiguous_bounds(), binned_cube.coord('latitude').contiguous_bounds()
    binned_data, _, _ = numpy.histogram2d(lons, lats, bins=[xs, ys])
    binned_cube.data = numpy.transpose(binned_data)
    return binned_cube

    
def storm_lats_lons(storms, years, months, basin, genesis=False, 
                 lysis=False, max_intensity=False):
    """ 
    Returns array of latitude and longitude values for all storms that 
    occurred within a desired year, month set and basin. 
    
    To get genesis, lysis or max intensity results set:
    Genesis plot: genesis=True
    Lysis plot: lysis=True
    Maximum intensity (location of max wind): 
    max_intensity=True
    
    """
    lats, lons = [], []
    count = 0
    for year in years:
        for storm in _storms_in_time_range(storms, year, months):
            if _storm_in_basin(storm, basin):
                if genesis:
                    #print 'getting genesis locations'
                    lats.extend([storm.obs_at_genesis().lat])
                    lons.extend([storm.obs_at_genesis().lon])
                elif lysis:
                    #print 'getting lysis locations'
                    lats.extend([storm.obs_at_lysis().lat])
                    lons.extend([storm.obs_at_lysis().lon])
                elif max_intensity:
                    #print 'getting max int locations'
                    lats.extend([storm.obs_at_vmax().lat])
                    lons.extend([storm.obs_at_vmax().lon])
                else:
                    #print 'getting whole storm track locations'
                    lats.extend([ob.lat for ob in storm.obs])
                    lons.extend([ob.lon for ob in storm.obs])
                count += 1
                
    # Normalise lon values into the range 0-360
    norm_lons = []
    for lon in lons:
        norm_lons.append((lon + 720) % 360)
    return lats, norm_lons, count


def get_projected_track(storm, map_proj):
    """ Returns track of storm as a linestring """
    lons, lats = list(zip(*[(ob.lon, ob.lat) for ob in storm.obs]))
    track = sgeom.LineString(list(zip(lons, lats)))
    projected_track = map_proj.project_geometry(track, ccrs.Geodetic())
    return projected_track

def _get_annual_vmax_storm_count_hemi(storms, years, months, basin):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = []
    count = 0
#    print years, months
    for storm in ts_model.example_code._storms_in_time_range(storms, years[0], months):
#        print 'found storm'
        if _storm_vmax_in_basin(storm, basin):
            count += 1
    storm_counts.append(count)
    return storm_counts    

def _get_annual_vmax_storm_ace_hemi(storms, years, months, basin):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_ace = []
    ace = 0
#    print years, months
    for storm in ts_model.example_code._storms_in_time_range(storms, years[0], months):
        #print 'found storm in ace, anywhere ', storm.ace_index()
        if _storm_vmax_in_basin(storm, basin):
            ace += storm.ace_index_no6hrcheck()
    storm_ace.append(ace)
    return storm_ace   

def get_annual_vmax_mean_count(storms, years, months, basin, nensemble):
    """ 
    Returns list of annual mean storm counts for a given
    set of years and months 
    
    
    """
    mean_annual_count = []
    for year in years:
        annual_count = _get_annual_vmax_storm_count_hemi(storms, [year,year], months, basin)
        for count in annual_count:
            mean_annual_count.append(float(count)/float(nensemble))
    return mean_annual_count

def get_annual_vmax_mean_ace(storms, years, months, basin, nensemble):
    """ 
    Returns list of annual mean ace for a given
    set of years and months 
    """
    mean_annual_ace = []
    for year in years:
        annual_ace = _get_annual_vmax_storm_ace_hemi(storms, [year,year], months, basin)
#        print 'year, annual ace ',annual_ace
        for ace in annual_ace:
            mean_annual_ace.append(float(ace)/float(nensemble))
    return mean_annual_ace

def annual_vmax_storm_counts(storms, years, months, basin, storm_types=['SS', 'TS', 'HU', 'MH']):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = []
    for year in years:
        count = 0
        for storm in ts_model.example_code._storms_in_time_range(storms, year, months):
            if (storm.max_storm_type() in storm_types) and _storm_vmax_in_basin(storm, basin):
                count += 1
        storm_counts.append(count)
    return storm_counts    

def annual_vmax_storm_ace(storms, years, months, basin, storm_types=['SS', 'TS', 'HU', 'MH']):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_ace = []
    for year in years:
        ace = 0
        for storm in ts_model.example_code._storms_in_time_range(storms, year, months):
            if (storm.max_storm_type() in storm_types) and _storm_vmax_in_basin(storm, basin):
                ace += storm.ace_index()
        storm_ace.append(ace)
    return storm_ace  

def storm_intensity(mslp):
    """ 
    Returns hurricane category based on Saffir-Simpson Hurricane 
    Wind Scale. Non-hurricanes return '--' 
    
    """
#    print 'mslp',mslp
    if (float(mslp) >= 994.):
        category=0 # Category 1
    elif (float(mslp) >= 980. and float(mslp) < 994.):
        category=1 # Category 2
    elif (float(mslp) >= 965. and float(mslp) < 980.):
        category=2 # Category 3  
    elif (float(mslp) >= 945. and float(mslp) < 965.):
        category=3 # Category 4       
    elif (float(mslp) >= 920. and float(mslp) < 945.):
        category=4 # Category 4       
    elif (float(mslp) >= 860. and float(mslp) < 920):
        category=5 # Category 5                   
    else:
        category=0
    return category
    
def storm_intensity_vmax(vmax):
    """ 
    Returns hurricane category based on Saffir-Simpson Hurricane 
    Wind Scale. Non-hurricanes return '--' 
    
    """
    if vmax >= 64 and vmax <= 82:
        category=1 # Category 1
    elif vmax >= 83 and vmax <= 95:
        category=2 # Category 2
    elif vmax >= 96 and vmax <= 112:
        category=3 # Category 3  
    elif vmax >= 113 and vmax <= 136:
        category=4 # Category 4       
    elif vmax >= 137:
        category=5 # Category 5                   
    else:
        category=0
    return category

