""" 
Plots the (optionally) tracks, genesis locations, lysis (dissipation)
locations and location of maximum intensity of model tropical
storms for a desired set of years and months.

Tracks are not coloured according to maximum intensity. 

In order to choose specific basins, you may need the python iris module (within th storm_assess/__init__.py code) in order to find if storms cross the basin.
This capability is currently commented out for simplicity
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os, sys, glob
import storm_assess
import storm_assess.functions
import storm_assess.track_cmor as track_cmor
import netCDF4, cftime
import numpy as np

# some example inputs for the plotting
years = np.arange(1950,2015)
months = [5,6,7,8,9,10,11]
months_nh = [5,6,7,8,9,10,11]
months_sh = [10,11,12,1,2,3,4,5]

# plot all basin storms
basin = None

experiment = 'highresSST-present'
YEARSTART = '1950'
YEAREND = '2014'
algorithm = 'TRACK'
yearsplot = list(range(YEARSTART, YEAREND+1))
    
institute = 'CMCC'
model = 'CMCC-CM2'
resol = 'HR4'
model_grid = 'gn'
member_id = 'r1i1p1f1'

institute = 'ECMWF'
model = 'ECMWF-IFS'
resol = 'LR'
model_grid = 'gr'
member_id = 'r1i1p1f1'

# directory containing the storm track files
base_dir = '/gws/nopw/j04/primavera1/model_derived_data/storm_tracking/'
CMOR_DIR = os.path.join(base_dir, algorithm, institute, model+'-'+resol, experiment, member_id, 'tropical', 'v4')

def storm_tracks(storms, years, months, basin, title, fig, ax, algorithm, hemi, genesis=False, lysis=False, max_intensity=False, warmcore = False, yoff=0.):
    """ 
    Plots storm tracks, genesis locations, lysis locations 
    or location of maximum intensity for all model tropical
    storms that form in a desired set of years, months and 
    ocean basin. Default plot is of storm tracks. 
    
    To get different plots set:
    Genesis plot: genesis=True
    Lysis plot: lysis=True
    Maximum intensity (location of max wind): 
    max_intensity=True
    
    Basin options:
    None: Whole globe
    'na': North Atlantic
    'ep': Eastern Pacific
    'wp': Western Pacific
    'ni': North Indian Ocean
    'si': Southwest Indian Ocean
    'au': Australian Region
    'sp': South Pacific
    
    Note: months [1,2,3,4,5,6] will obtain storms 
    that *formed* within the time period 1 Jan to
    30 June inclusive. 
    
    Setting months [11,12,1,2,3,4] and years [1996]
    will return storms that formed between 1 Nov
    1996 and 30 April 1997 inclusive, assuming those 
    data are available in the input file.
    
    Setting the basin will return any storms which
    *passed through* (not just formed) in the 
    designated basin region. 
    
    
    """   
    count = 0
    for year in years:
        for storm in storm_assess.functions._storms_in_time_range(storms, year, months):
            if genesis:
                variable = 'Genesis'
                ax.plot(storm.obs_at_genesis().lon, storm.obs_at_genesis().lat,
                         'bo', markersize=3, transform=ccrs.Geodetic())

            elif lysis:
                variable = 'Lysis'
                ax.plot(storm.obs_at_lysis().lon, storm.obs_at_lysis().lat,
                         'go', markersize=3, transform=ccrs.Geodetic())

            elif max_intensity:
                variable = 'Maximum Intensity'
                ax.plot(storm.obs_at_vmax().lon, storm.obs_at_vmax().lat,
                         'ro', markersize=3, transform=ccrs.Geodetic())

            else:
                variable = 'Tracks'   
                ax.plot([ob.lon for ob in storm.obs], [ob.lat for ob in storm.obs],
                     linewidth=1.2, transform=ccrs.Geodetic())
            count += 1

    if count != 0:
        fig.gca().coastlines() 
        title1 = 'Model Tropical Storm %s\n %s (%s - %s) using %s \n %s' % \
                  (variable, storm_assess.functions.BASIN_NAME.get(basin), str(years[0]), str(years[-1]), algorithm, title)
        print(title1)
        ax.set_title(title1, fontsize=12)
        s = ('Total tropical storms for %s: %s' % (hemi, count))
        fig.text(0.02, 0.09-yoff, s, ha='left', va='center', transform=plt.gca().transAxes)
        #print s   
    else:
        print('No storms found')

def read_storms(dir_in, hemi, run, yearstart, yearend):
    '''
    For a given file name format, and directory
    Read the storms from the files
    '''

    fname_nc = 'TC-{}_{}_{}-{}_{}_{}_{}_{}0101-{}1231.nc*'
   
    search_nh = fname_nc.format('NH', run['algorithm'], run['model'], run['resol'], experiment, member_id, run['grid'], yearstart, yearend)
    fname_nh = glob.glob(os.path.join(dir_in, search_nh))
    search_sh = fname_nc.format('SH', run['algorithm'], run['model'], run['resol'], experiment, member_id, run['grid'], yearstart, yearend)
    fname_sh = glob.glob(os.path.join(dir_in, search_sh))

    if hemi == 'nh':
        months = months_nh
        if len(fname_nh) > 0:
            fname = fname_nh[0]
        else:
            raise Exception('NH filename not found '+search_nh)
    else:
        months= months_sh
        if len(fname_sh) > 0:
            fname = fname_sh[0]
        else:
            raise Exception('NH filename not found '+search_sh)

    # derive path to data and read the netcdf file
    #path = os.path.join(dir_in, fname)
    path = fname
    print('path ',path)
    storms = list(track_cmor.load_cmor(path))
    print('no storms ',len(storms))

    # check the metadata to discover which algorithm this is, and hence
    # what feature variable is tracked
    with netCDF4.Dataset(path, 'r') as nc:
        track_algorithm = nc.getncattr('algorithm')
        if track_algorithm == 'TRACK':
            track_extra = nc.getncattr('algorithm_extra')
            if track_extra == 'T63avg':
                feature_variable = 'vortmean_T63'
            else:
                feature_variable = 'rv850_T42'
        elif track_algorithm == 'TempestExtremes':
            feature_variable = 'psl'
        else:
            raise Exception('Unrecognised algorithm in netcdf file '+path)

    return storms, feature_variable, track_extra

def work(runid_info, data_dir, yearstart, yearend, yearsplot):

    fig = plt.figure(figsize=(9,6), dpi=100)
    ax = storm_assess.functions.load_map(basin=basin)

    yoff = 0
    for hemi in ['nh', 'sh']:
        storms, feature_variable, track_extra = read_storms(data_dir, hemi, runid_info, yearstart, yearend)

        #for storm in storms:
        #    print storm.genesis_date()

        if hemi == 'sh': yoff = 0.05

        title = runid_info['model']+'-'+runid_info['resol']+', '+str(yearstart)+'-'+str(yearend)
        storm_tracks(storms, yearsplot, months, basin, title, fig, ax, algorithm, hemi, yoff=yoff, max_intensity = False, genesis = False)

    current_dir = os.getcwd()
    figname = os.path.join(current_dir, runid_info['model']+algorithm+'.png')
    plt.savefig(figname)

    print(years)
    plt.show()
    
def test_360day():
    '''
    Set up some test data values
    The test data is in the test_data subdirectory
    '''
    institute = 'MOHC'
    model = 'HadGEM3-GC31'
    resol = 'LM'
    model_grid = 'gn'
    member_id = 'r1i1p1f1'
    experiment = 'highresSST-present'
    algorithm = 'TRACK'
    yearstart = '2014'
    yearend = '2014'
    yearsplot = list(range(2014, 2015))
    runid_info = {'model': model, 'resol': resol, 'grid': model_grid, 'algorithm': algorithm}

    # the test data is in this subdirectory 
    dir_test = os.path.join(os.getcwd(), 'test_data')

    # call the read and plot subroutine
    work(runid_info, dir_test, yearstart, yearend, yearsplot)

def test_gregorian():
    '''
    Set up some test data values
    The test data is in the test_data subdirectory
    '''
    institute = 'ECMWF'
    model = 'ECMWF-IFS'
    resol = 'LR'
    model_grid = 'gr'
    member_id = 'r1i1p1f1'
    experiment = 'highresSST-present'
    algorithm = 'TRACK'
    yearstart = '1950'
    yearend = '2014'
    yearsplot = list(range(1980, 1981))
    runid_info = {'model': model, 'resol': resol, 'grid': model_grid, 'algorithm': algorithm}

    # the test data is in this subdirectory 
    dir_test = os.path.join(os.getcwd(), 'test_data')

    # call the read and plot subroutine
    work(runid_info, dir_test, yearstart, yearend, yearsplot)

def test_noleap():
    '''
    Set up some test data values
    The test data is in the test_data subdirectory
    '''
    institute = 'CMCC'
    model = 'CMCC-CM2'
    resol = 'HR4'
    model_grid = 'gn'
    member_id = 'r1i1p1f1'
    experiment = 'highresSST-present'
    algorithm = 'TRACK'
    yearstart = '2014'
    yearend = '2014'
    yearsplot = list(range(2014, 2015))
    runid_info = {'model': model, 'resol': resol, 'grid': model_grid, 'algorithm': algorithm}

    # the test data is in this subdirectory 
    dir_test = os.path.join(os.getcwd(), 'test_data')

    # call the read and plot subroutine
    work(runid_info, dir_test, yearstart, yearend, yearsplot)

if __name__ == '__main__':

    #test_360day()
    test_gregorian()
    #test_noleap()

    #runid_info = {'model': model, 'resol': resol, 'grid': model_grid, 'algorithm': algorithm}
    #work(runid_info, CMOR_DIR, YEARSTART, YEAREND)
