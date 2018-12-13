""" 
Load function to read in Reading Universities TRACK output. 
Written for the CMORised netcdf TRACK files
Should work for both an individual month and all months 


"""
import os.path, sys
import collections
import datetime
import netcdftime
import netCDF4

""" Store model storm observations as a named tuple (this enables access to data by its 
name instead of position index) """
Observation = collections.namedtuple('Observation', ['date', 'lat', 'lon', 'vort', 'vmax', 
                                                     'mslp', 't63_850', 't63_700', 't63_600', 't63_500', 't63_250', 't63_850_250_diff', 'extras'])

class Observation(Observation):
    """  Represents a single observation of a model tropical storm. """
    
    def six_hourly_timestep(self):
        """ Returns True if a storm record is taken at 00, 06, 12 or 18Z only """
        return self.date.hour in (0, 6, 12, 18) and self.date.minute == 0 and self.date.second == 0
    
    def add_to_axes(self, ax):
        """ Instructions on how to plot a model tropical storm observation """
        ax.plot(self.lon, self.lat)
        

class Storm(object):
    def __init__(self, snbr, obs, extras=None):
        """ Stores information about the model storm (such as its storm number) and corresponding 
        observations. Any additional information for a storm should be stored in extras as a dictionary. """
        self.snbr = snbr
        self.obs = obs
        
        # Create an empty dictionary if extras is undefined
        if extras is None:
            extras = {}
        self.extras = extras 

    @property
    def vmax(self):
        """ The maximum wind speed attained by the storm during its lifetime """
        return max(ob.vmax for ob in self.obs)
    
    @property
    def mslp_min(self):
        """ The minimum central pressure reached by the storm during its lifetime (set to 
        -999 if no records are available) """
        mslps = [ob.mslp for ob in self.obs if ob.mslp != 1e12]  
        if not mslps:
            mslps = [-999]
        return min(mslps)
    
    @property
    def vort_max(self):
        """ The maximum 850 hPa relative vorticity attained by the storm during its lifetime """
        return max(ob.vort for ob in self.obs)
    
    @property
    def t63_850_max(self):
        """ The maximum T63 850 vorticity attained by the storm during its lifetime """
        return max(ob.t63_850 for ob in self.obs)
    
    @property
    def t63_250_max(self):
        """ The maximum T63 250 vorticity attained by the storm during its lifetime """
        return max(ob.t63_250 for ob in self.obs)
    
    @property
    def t63_850_250_diff_max(self):
        """ The maximum T63 850 - 250 vorticity attained by the storm during its lifetime """
        return max(ob.t63_850_250_diff for ob in self.obs)
    
    def __len__(self):
        """ The total number of observations for the storm """
        return len(self.obs)
    
    def nrecords(self):
        """ The total number of records/observations for the storm """
        return len(self)
    
    def number_in_season(self):
        """ Returns storm number of the storm (the number of the storm for that year and ensemble
        member number) """
        return self.snbr    
    
    def lifetime(self):
        """ The total length of time that the storm was active. This uses all observation
        points, no maximum wind speed threshold has been set """
        return max(ob.date for ob in self.obs)-min(ob.date for ob in self.obs)
        
    def genesis_date(self):
        """ The first observation date that a storm becomes active """
        #return min(ob.date for ob in self.obs)
        return self.obs_at_genesis().date
    
    def lysis_date(self):
        """ The final date that a storm was active """
        #return max(ob.date for ob in self.obs)
        return self.obs_at_lysis().date
    
    def ace_index(self):
        """ The accumulated cyclone energy index for the storm. Calculated as the square of the
        storms maximum wind speed every 6 hours (0, 6, 12, 18Z) throughout its lifetime. Observations
        of the storm taken in between these records are not currently used. Returns value rounded to
        2 decimal places. Wind speed units: knots """
        ace_index = 0
        for ob in self.obs:
            if ob.six_hourly_timestep():
                ace_index += numpy.square(ob.extras['vmax_kts'])/10000.
        return round(ace_index, 2)
    
    def obs_at_vmax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return max(self.obs, key=lambda ob: ob.vmax)
    
    def obs_at_vortmax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return max(self.obs, key=lambda ob: ob.vort)
    
#    def obs_at_coremax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
#	t63_diff = [ob.t63_diff for ob in self.obs]
#	peak = find_local_maximum(t63_diff)
        #return max(self.obs, key=lambda ob: ob.t63_diff)
#        return [ob for ob in self.obs][peak]
    
    def obs_at_mslpmin(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return min(self.obs, key=lambda ob: ob.mslp)
    
    def obs_at_genesis(self):
        """Returns the Observation instance for the first date that a storm becomes active """       
        for ob in self.obs:
            return ob
        else:
            raise ValueError('model storm was never born :-(')

    def obs_at_lysis(self):
        """Returns the Observation instance for the last date that a storm was active """    
        return [ob for ob in self.obs][-1]
    

def load_cmor(fh):
    '''
    Load cmor netcdf format of track files
    The input file should contain (at least):
       Dimensions:
          ntracks: number of storms in file
          record: total number of time points
          plev: number of pressure levels (with associated data)
       Variables:
          TRACK_ID(tracks): Storm number
          FIRST_PT(tracks): Index to first point in each storm
          NUM_PTS(tracks): The number of points in this storm
          index(record): storm track sequence number (index to this storm in the whole record diemsion)
          vortmean_T63(record): feature tracked variable
          lon(record): longitude of feature tracked variable
          lat(record): latitude of feature tracked variable
    '''
    # final variables needed
    # storm number snbr
    # date, lat, long, vort, vmax, mslp, T63[nlev], vmax_kts, w10m

    scaling_ms_knots = 1.944
    print 'fh type', type(fh)
    fh_type = str(type(fh))
    if 'str' in fh_type:
        fh = [fh]
    
    # for each file in the file handle            
    for fname in fh:
        if not os.path.exists(fname):
            raise Exception('Input file does not exist '+fname)
        else:
            print 'fname ',fname
            with netCDF4.Dataset(fname, 'r') as nc:
                track_algorithm = nc.getncattr('algorithm')
                if track_algorithm == 'TRACK':
                    if nc.getncattr('algorithm_extra') == 'T63avg':
                        vort_variable = 'vortmean_T63'
                    else:
                        vort_variable = 'rv850_T42'
                elif track_algorithm == 'TempestExtremes':
                    vort_variable = 'rv850'

                # number of storms in the file
                ntracks = int(nc.dimensions['tracks'].size)
                plev = int(nc.dimensions['plev'].size)
        # Loop through each storm, and create a class object containing the storm properties
                psl = nc.variables['psl']
                if psl.units == 'Pa':
                    psl_scaling = 1.0 / 1000.0
                else:
                    psl_scaling = 1.0

                # read the time variable and convert to a more useful format
                time_var = nc.variables['time']
                dtime = netCDF4.num2date(time_var[:],time_var.units, calendar = time_var.calendar)

                first_pts = nc.variables['FIRST_PT']
                storm_lengths = nc.variables['NUM_PTS']
                indices = nc.variables['index']
                lats = nc.variables['lat']
                lons = nc.variables['lon']
                vorts = nc.variables[vort_variable]
                psls = nc.variables['psl']
                sfcWinds = nc.variables['sfcWind']
                vmaxs = nc.variables['ws925']
                
                # number of pressure levels for variables
                if plev >= 5:
                    rv850_T63 = nc.variables['rv850_T63']
                    rv700_T63 = nc.variables['rv700_T63']
                    rv600_T63 = nc.variables['rv600_T63']
                    rv500_T63 = nc.variables['rv500_T63']
                    rv250_T63 = nc.variables['rv250_T63']
                    rv_diff_850_250 = (rv850_T63[:] - rv250_T63[:])
                elif plev < 5:
                    rv850_T63 = nc.variables['rv850_T63']
                    rv500_T63 = nc.variables['rv500_T63']
                    rv250_T63 = nc.variables['rv250_T63']
                    rv_diff_850_250 = (rv850_T63[:] - rv250_T63[:])

                for storm_no in range(ntracks):
                #for storm_no in range(2):
                    storm_obs = []
                    tcid = storm_no
                    #print 'storm_no ',storm_no
                    first_pt = first_pts[storm_no]
                    storm_length = storm_lengths[storm_no]
                    record_no = storm_length
                    index = indices[first_pt:first_pt+storm_length]
                    #print storm_no, index
                    for ip in index:
                        i = ip+first_pt
                        date = dtime[i]
                        #print storm_no, ip, i, first_pt, date

                        # lat, lon are the latitude and longitude associated with the tracked feature
                        lat = lats[i]
                        lon = lons[i]
                        vort = vorts[i]
                        psl = psls[i] * psl_scaling
                        sfcWind = sfcWinds[i]
                        vmax = vmaxs[i]
                        vmax_kts = vmax * scaling_ms_knots

                        rv850_t63_this = rv850_T63[i]
                        rv500_t63_this = rv500_T63[i]
                        rv250_t63_this = rv250_T63[i]
                        rv_diff_850_250_this = rv_diff_850_250[i]
                        if plev >= 5:
                            rv700_t63_this = rv700_T63[i]
                            rv600_t63_this = rv600_T63[i]
                        else:
                            rv700_t63_this = -99.9
                            rv600_t63_this = -99.9
            
                        storm_obs.append(Observation(date, lat, lon, vort, vmax, psl, rv850_t63_this, rv700_t63_this, rv600_t63_this, rv500_t63_this, rv250_t63_this, rv_diff_850_250_this, extras={'vmax_kts':vmax_kts, 'w10m':sfcWind}))
                    #print 'yielding ',storm_obs

            # Yield storm
                    yield Storm(tcid, storm_obs, extras={})

if __name__ == '__main__':
    fname = os.path.join(SAMPLE_TRACK_DATA)
    print 'Loading TRACK data from file:' , fname    
    storms = list(load(fname, ex_cols=3, calendar='netcdftime'))
    print 'Number of model storms: ', len(storms)
    
    # Print storm details:
    for storm in storms: 
        #print storm.snbr, storm.genesis_date()
        for ob in storm.obs:
            print ob.date, ob.lon, ob.lat, ob.vmax, ob.extras['vmax_kts'], ob.mslp, ob.vort
    print 'Number of model storms: ', len(storms)
    
