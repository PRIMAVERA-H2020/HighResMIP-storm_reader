import numpy
import os.path
import collections

import datetime
import fnmatch

import cartopy.crs as ccrs
import shapely.geometry as sgeom
import scipy.signal as signal

def find_local_maximum(data):
    max_peakind = signal.find_peaks_cwt(data, numpy.arange(1,100))
    print 'max_peak ',max_peakind
    print 'data ',data
    for peak in max_peakind:
        local_max = data[peak]
	print 'peak ',peak
	print 'local_max ',local_max, data[max([peak-4, 0]):peak]
	min_index = max([peak-4, 0])
	max_index = min([peak+5,len(data)])
        if max(data[min_index:peak]) <= local_max and \
                   local_max >=  max(data[peak+1:max_index]):
            return peak

def _basin_polygon(basin, project=True):
    # Lat/lon locations of tracking regions for each ocean basin
    # note that these are not necessarily those defined by WMO
    TRACKING_REGION = {'na': ([-75, -20, -20, -80, -80, -100, -100, -75, -75], [0, 0, 60, 60, 40, 40, 20, 6, 0]),
                       #'ep': ([-140, -75, -75, -100, -100, -120, -120, -140, -140], [0, 0, 6, 20, 40, 40, 60, 60, 0]),
                       'ep': ([-140, -75, -75, -100, -100, -140, -140], [0, 0, 6, 20, 30, 30, 0]),
                       'wp': ([-260, -180, -180, -260, -260], [0, 0, 60, 60, 0]),
                       #'wp': ([-260, -180, -180, -260, -260], [0, 0, 25, 25, 0]),
                       'ni': ([-320, -260, -260, -320, -320], [0, 0, 30, 30, 0]),
                       'si': ([-330, -270, -270, -330, -330], [-40, -40, 0, 0, -40]),
                       'au': ([-270, -200, -200, -270, -270], [-40, -40, 0, 0, -40]),
                       'sp': ([-200, -120, -120, -200, -200], [-40, -40, 0, 0, -40]),
                       'sa': ([-90, 0, 0, -90, -90], [-40, -40, 0, 0, -40]),
                       #'mdr': ([-80, -20, -20, -80, -80], [10, 10, 20, 20, 10])
                       'mdr': ([-60, -20, -20, -60, -60], [10, 10, 20, 20, 10])
                       }
    
    rbox = sgeom.Polygon(zip(*TRACKING_REGION.get(basin)))
    if project: 
        rbox = ccrs.PlateCarree().project_geometry(rbox, ccrs.PlateCarree())
    return rbox
    
    
def _storm_in_basin(storm, basin):
    """ Returns True if a storm track intersects a defined ocean basin """
    rbox = _basin_polygon(basin)   
    lons, lats = zip(*[(ob.lon, ob.lat) for ob in storm.obs])
    track = sgeom.LineString(zip(lons, lats))       
    projected_track = ccrs.PlateCarree().project_geometry(track, ccrs.Geodetic())
    if rbox.intersects(projected_track):
        return True
    return False

def _storm_vmax_in_basin(storm, basin):
    """ Returns True if the maximum intensity of the storm occurred
    in desired ocean basin. """
    rbox = _basin_polygon(basin)  
    xy = ccrs.PlateCarree().transform_point(storm.obs_at_vmax().lon, storm.obs_at_vmax().lat, ccrs.Geodetic())
    point = sgeom.Point(xy[0], xy[1])
    if point.within(rbox):
        return True
    return False

def _storm_genesis_in_basin(storm, basin):
    """ Returns True if the maximum intensity of the storm occurred
    in desired ocean basin. """
    rbox = _basin_polygon(basin)  
    xy = ccrs.PlateCarree().transform_point(storm.obs_at_genesis().lon, storm.obs_at_genesis().lat, ccrs.Geodetic())
    point = sgeom.Point(xy[0], xy[1])
    if point.within(rbox):
        return True
    return False
    
def _storms_in_basin_year_member_forecast(storms, basin, years, members, fcst_dates):
    """ 
    A generator which yields storms that occur within a desired ocean basin 
    with a particular start date, ensemble member number and forecast date 
    
    """
    for storm in storms:
        if (storm.genesis_date().year in years) and \
            (storm.extras['member'] in members) and \
            (storm.extras['fcst_start_date'] in fcst_dates) and \
            _storm_in_basin(storm, basin):
            yield storm  

def _storms_in_year_member_forecast(storms, years, members, fcst_dates):
    for storm in storms:
        if (storm.genesis_date().year in years) and \
            (storm.extras['member'] in members) and \
            (storm.extras['fcst_start_date'] in fcst_dates):
            yield storm
             
             
#def _get_time_range(year, months):
#    """ 
#    Creates a start and end date (a datetime.date timestamp) for a 
#    given year and a list of months. If the list of months overlaps into 
#    the following year (for example [11,12,1,2,3,4]) then the end date 
#    adds 1 to the original year 
#    
#    """
#    start_date = datetime.datetime(year, months[0], 1)
#    end_year = year
#    end_month = months[-1]+1
#    if months[-1]+1 < months[0] or months[-1]+1 == 13 or len(months) >= 12:
#        end_year = year+1
#    if months[-1]+1 == 13:
#        end_month = 1
#    end_date = datetime.datetime(end_year, end_month, 1)
#    return start_date, end_date
#                
        
#def _storms_in_time_range(storms, year, months):
#    """Returns a generator of storms that formed during the desired time period """
#    start_date, end_date = _get_time_range(year, months)
#    for storm in storms:        
#        if (storm.genesis_date() >= start_date) and (storm.genesis_date() < end_date):
#            yield storm
            
            
def _storms_in_basin_year_month_member_forecast(storms, basin, year, months, members, fcst_dates):
    """ 
    A generator which yields storms that occurred within a desired ocean basin 
    with a particular start date, start month, ensemble member number and forecast date 
    
    """
    for storm in _storms_in_time_range(storms, year, months):
        if (storm.extras['member'] in members) and \
           (storm.extras['fcst_start_date'] in fcst_dates) and \
            _storm_in_basin(storm, basin):
            yield storm
            
          
#def _tracking_file_exists(start_date, year, member, hemisphere):
#    """ 
#    Searches for a filename with a given forecast start date, year
#    and ensemble member number. Returns True if file is found.
#    
#    """
#    for root, dirs, files in os.walk(TRACK_DIR):
#        for file in files:
#            fname = os.path.join(root, file)
#            if os.path.isfile(fname):
#                if fnmatch.fnmatch(fname, 
#                                   '%sff_trs.vor_fullgrid_wind_mslp_L5.new.%s_%s_%s.%s.date' % 
#                                   (TRACK_DIR, str(year), start_date, member, hemisphere)
#                                   ):
#                    return True
             
def _tracking_file_exists(fcst_date, year, member, hemisphere, 
                          model='glosea5', file_type='tropical'):
    """ 
    Searches for a filename with a given forecast start date, year
    and ensemble member number. Returns True if file is found.
    
    """
    tracking_dir = TRACK_DIR + model
    
    for root, dirs, files in os.walk(tracking_dir):
        for file in files:
            fname = os.path.join(root, file)
            if os.path.isfile(fname):
                if fnmatch.fnmatch(fname, 
                    '%s/ff_trs.vor_fullgrid_wind_mslp_L5.new.%s.%s_%s_%s.%s.date' % 
                    (tracking_dir, file_type, str(year), fcst_date, 
                     member, hemisphere)):
                    return True
                
                
def ensemble_count(years, fcst_dates, members, hemisphere, file_type='tropical'):
    """ Returns number of individual tracking files available for a given 
    set of forecast dates, ensemble members and years """
    count  = 0
    for year in years:
        for fcst_date in fcst_dates:
            for member in members:
                if _tracking_file_exists(fcst_date, year, member, 
                                         hemisphere, file_type=file_type):
                    count += 1
    return count

