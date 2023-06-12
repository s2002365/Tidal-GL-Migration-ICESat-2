##############################################################################
############ FUNCTIONS TO READ AND PROCESS ATL06 DATA FOR RTLA ###############
############ Author: Bryony Freer                              ###############
##############################################################################

# We acknowledge the code written in ICESat-2 hackweek tutorials, for which some of this code is based from


import glob
import os
import sys
import re
import pyproj
import datetime as dt
import numpy as np
from numpy import arange, array, nonzero, in1d, abs, linspace, sqrt, pi
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import h5py
import pyproj
import gdal
import osr
from scipy.interpolate import interp1d
from scipy import interpolate

try:
    import pointCollection as pc
except Exception:
    import pointCollection as pc

import ipywidgets as widgets
from IPython.display import HTML

from pyTMD.read_tide_model import extract_tidal_constants
from pyTMD.read_netcdf_model import extract_netcdf_constants
from pyTMD.read_GOT_model import extract_GOT_constants
from pyTMD.read_FES_model import extract_FES_constants
from pyTMD.compute_tide_corrections import compute_tide_corrections

def create_color_dict():
    return {3:'#1f77b4', 4: '#ff7f0e', 5: '#2ca02c', 6: '#d62728', 7: '#9467bd', 8: '#8c564b', 9: '#e377c2', 10: '#7f7f7f', 11:'#bcbd22', 12: '#17becf', 13:'#6b2472', 14:'#55FF55', 15:'#f53b92', 16: '#2D2B83', 17:'#FBFF00'}

#Define function to calculate distance between coordinates
def distance(lat1,lon1,lat2,lon2):
    '''
    Calculate distance between a set point and a lat lon pair from the data
    lat1,lon1 = location 1 
    lat2,lon2 = location 2
    '''
    
    # approximate radius of earth in km
    R = 6373.0

    lat1r = radians(lat1)
    lon1r = radians(lon1)
    lat2r = radians(lat2)
    lon2r = radians(lon2)

    dlon = lon2r - lon1r
    dlat = lat2r - lat1r

    a = sin(dlat / 2)**2 + cos(lat1r) * cos(lat2r) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))

    distance = R * c
    return distance

def atl06_to_dict(filename, beam, field_dict=None, index=None, epsg=None):
    """
        Read selected datasets from an ATL06 file

        Input arguments:
            filename: ATl06 file to read
            beam: a string specifying which beam is to be read (ex: gt1l, gt1r, gt2l, etc)
            field_dict: A dictinary describing the fields to be read
                    keys give the group names to be read, 
                    entries are lists of datasets within the groups
            index: which entries in each field to read
            epsg: an EPSG code specifying a projection (see www.epsg.org).  Good choices are:
                for Greenland, 3413 (polar stereographic projection, with Greenland along the Y axis)
                for Antarctica, 3031 (polar stereographic projection, centered on the Pouth Pole - same as MOA)
        Output argument:
            D6: dictionary containing ATL06 data.  Each dataset in 
                dataset_dict has its own entry in D6.  Each dataset 
                in D6 contains a numpy array containing the 
                data. Keys: 
                - 'latitude': array
                - 'longitude': array
                - 'h_li': array
                - 'atl06_quality_summary': array
                - 'x_atc'
                - 'y_atc': array
                - 'dh_fit_dx': array
                - 'dh_fit_dy': array
                - 'x': array
                - 'y': array
                - 'rgt': integer (reference ground track number)
                - 'cycle': integer
                - 'beam': string of beam e.g. '/gt2l'

    """
    if field_dict is None:
        field_dict={None:['latitude','longitude','h_li', 'atl06_quality_summary'],\
                    'ground_track':['x_atc','y_atc'],\
                    'fit_statistics':['dh_fit_dx', 'dh_fit_dy']}
    
    # Create empty dictionary to fill with subsequent information in correct structure.
    D={}
    
    # Following line tears apart the ATL06 naming convention and tells you what the different parts of the filename mea (date, rgt, cycle, region, release, version). Stores all in dictionaries.
    file_re=re.compile('ATL06_(?P<date>\d+)_(?P<rgt>\d\d\d\d)(?P<cycle>\d\d)(?P<region>\d\d)_(?P<release>\d\d\d)_(?P<version>\d\d).h5')
    
    # Use H5PY module to open the file 
    with h5py.File(filename,'r') as h5f:
        for key in field_dict:
            for ds in field_dict[key]:
                #For each beam in the field dictionary, pull out the variables in land ice segments, and stores them by their field name in the dictionary (as numpy arrays).
                if key is not None:
                    ds_name=beam+'/land_ice_segments/'+key+'/'+ds
                else:
                    ds_name=beam+'/land_ice_segments/'+ds
                if index is not None:
                    D[ds]=np.array(h5f[ds_name][index])
                else:
                    D[ds]=np.array(h5f[ds_name])
                if '_FillValue' in h5f[ds_name].attrs:
                    bad_vals=D[ds]==h5f[ds_name].attrs['_FillValue']
                    D[ds]=D[ds].astype(float)
                    D[ds][bad_vals]=np.NaN
    
    # Use PyProj to convert the lat/longs into x and y coordinates. 
    if epsg is not None:
        #xy=np.array(pyproj.proj.Proj(epsg)(D['longitude'], D['latitude'])) #original code that fails
        xy=np.array(pyproj.Proj("epsg:3031", preserve_units=False)(D['longitude'], D['latitude'])) #REMOVED EXTRA .proj, removed init=

        D['x']=xy[0,:].reshape(D['latitude'].shape)
        D['y']=xy[1,:].reshape(D['latitude'].shape)
    
    # Store cycle, rgt and beam information in the same dictionary
    temp=file_re.search(filename)
    D['rgt']=int(temp['rgt'])
    D['cycle']=int(temp['cycle'])
    D['beam']=beam
    return D

###############################################
## APPLY CROSS-TRACK SLOPE CORRECTION        ##
## And plot to see what the data looks like  ##
## and select cycles for RTLA                ##
###############################################

def min_seg_difference(D):
    """
    seg_difference_filter: Use elevations and slopes to find bad ATL06 segments
    
    Inputs: 
        D: a granule of ATL06 data, in dictionary format.  Must have entries:
            x_atc, h_li, dh_fit_dx
        
    Returns:
        delta_h_seg: the minimum absolute difference between each segment's endpoints and those of its two neighbors
    """
    h_ep=np.zeros([2, D['h_li'].size])+np.NaN
    h_ep[0, :]=D['h_li']-D['dh_fit_dx']*20
    h_ep[1, :]=D['h_li']+D['dh_fit_dx']*20
    delta_h_seg_pre=np.zeros_like(D['h_li'])
    delta_h_seg=np.zeros_like(D['h_li'])
    delta_h_seg[1:]=np.abs(D['h_li'][1:]-h_ep[1, :-1])
    delta_h_seg[:-1]=np.minimum(delta_h_seg[:-1], np.abs(D['h_li'][:-1]-h_ep[0, 1:]))
    return delta_h_seg

def PlotCorrectedRepeatTracks(data_root, rgt, rpt, color_dict, cross_track_corr=True, plot=True):
    
    '''Function to plot along-track elevation profile for each cycle of a single track. 
    Elevation values corrected for cross-track slope.'''
    
    #Create dictionary for L and R ground track of beam pair
    D_l={}
    D_r={}
    
    for cycle in ['03','04','05','06','07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17']:
        for filename in glob.glob(os.path.join(data_root, f'*processed_ATL06_*_{rgt}{cycle}*.h5')): 
            try:
                # read the left-beam data
                D_l[filename]=atl06_to_dict(filename,f'/gt{rpt}l', index=None, epsg=3031)
                # read the right-beam data
                D_r[filename]=atl06_to_dict(filename,f'/gt{rpt}r', index=None, epsg=3031)
            except Exception as e:
                print(f'filename={filename}, exception={e}')

    if plot == True:
        #Set up figure to plot the repeat track elevations
        fig=plt.figure(figsize=(8,8))
        ax=fig.add_subplot(111, ylabel='Elevation (m), Corrected \nfor Cross-Track Slope', xlabel='Latitude')
    
    for filename, Dl in D_l.items():
        # Also read in data for R beam 
        Dr=D_r[filename]

        # Set across-track correction parameter depending on RGT
        if rpt == 1:
            atc_corr_l=(Dl['y_atc']-3350)*Dl['dh_fit_dy']
            atc_corr_r=(Dr['y_atc']-3350)*Dr['dh_fit_dy']

        elif rpt == 2: 
            atc_corr_l=Dl['y_atc']*Dl['dh_fit_dy']
            atc_corr_r=Dr['y_atc']*Dr['dh_fit_dy']

        elif rpt == 3: 
            atc_corr_l=(Dl['y_atc']+3350)*Dl['dh_fit_dy']
            atc_corr_r=(Dr['y_atc']+3350)*Dr['dh_fit_dy']
        
        # Calculate minimum segment difference (want to remove all with values > 1, following hackweek tutorials, Arendt et al., 2019) 
        delta_h_seg_l = min_seg_difference(Dl)    
        delta_h_seg_r = min_seg_difference(Dr)
            
        # L beam: apply cross-track correction and quality filtering  --> add corrected elevations to dictionary      
        #correction 1 (apply to h_li)
        good=(Dl['atl06_quality_summary']==0) & (~(np.isnan(Dl['dh_fit_dy']))) & (delta_h_seg_l < 1)
        Dl['h_corr'] = Dl['h_li'][good]-atc_corr_l[good]
        Dl['lat_corr'] = Dl['latitude'][good] #adds latitude field with only values corresponding to high-qual, h_corr values
        Dl['lon_corr'] = Dl['longitude'][good]
        Dl['good_idx'] = good

        # R beam: (as above)
        good=(Dr['atl06_quality_summary']==0) & (~(np.isnan(Dr['dh_fit_dy']))) & (delta_h_seg_r < 1) 
        Dr['h_corr'] = Dr['h_li'][good]-atc_corr_r[good]
        Dr['lat_corr'] = Dr['latitude'][good]
        Dr['lon_corr'] = Dr['longitude'][good]
        Dr['good_idx'] = good
        
        if cross_track_corr == False:
            #L beam
            good=(Dl['atl06_quality_summary']==0) & (~(np.isnan(Dl['dh_fit_dy']))) & (delta_h_seg_l < 1)
            Dl['h_corr'] = Dl['h_li'][good]
            Dl['lat_corr'] = Dl['latitude'][good] #adds latitude field with only values corresponding to high-qual, h_corr values
            Dl['lon_corr'] = Dl['longitude'][good]
            Dl['good_idx'] = good
            
            #R beam
            good=(Dr['atl06_quality_summary']==0) & (~(np.isnan(Dr['dh_fit_dy']))) & (delta_h_seg_r < 1) 
            Dr['h_corr'] = Dr['h_li'][good]
            Dr['lat_corr'] = Dr['latitude'][good]
            Dr['lon_corr'] = Dr['longitude'][good]
            Dr['good_idx'] = good
            

        # Plot cross-track slope corrected data for L and R beams (filtered to include only high quality data)
        if plot == True:
            ax.plot(Dl['lat_corr'], Dl['h_corr'], '.', markersize=3,label=f"{Dl['cycle']} (L)", color=color_dict[int(Dl['cycle'])])
            ax.plot(Dr['lat_corr'], Dr['h_corr'], '+', markersize=4, label=f"{Dr['cycle']} (R)", color=color_dict[int(Dl['cycle'])])
        
    if plot == True:
        fig.suptitle(f'ATL06 Repeat Track Elevations: RGT {rgt}, Beam GT{rpt}\n (filtered using atl06_quality_summary parameter)\n', weight='bold')
        ax.title.set_text('With Cross-Track Slope Correction')
        if cross_track_corr == False:
            ax.title.set_text('Without Cross-Track Slope Correction')
        ax.legend(title='Cycle (Beam):', loc='upper left', bbox_to_anchor=(0, -0.2), ncol=2)
        plt.tight_layout();
        
    else: 
        print('Plotting turned off. Corrections applied to data.')
    
    return D_l, D_r



###############################################
## EXTRACT DATA INTO DICT OF REPEATS ##########
###############################################

def ExtractData(D, cycles, rgt, rpt, lr, atl06_dir,  tides=True):
    
    '''Function to extract the necessary data from D_l and D_r dictionaries as well as from along-track tides
    
    Inputs:
    - D: dictionary of values (either D_l or D_r, output of part 4)
    - cycles: list of cycles to include
    - rgt: reference ground track number as a four digit string e.g. "0034"
    - rpt: reference pair track number (1, 2 or 3)
    - lr: string "l" or "r" depending on left or right beam
    - atl06_dir: file path to directory containing processed ATL06 and CATS2008 along-track tide .h5 files
    - tides: set to True if want to extract tide data into the same dictionary.
    
    Outputs: 
    - dict_repeats: dictionary of values per cycle 
    '''
    
    #Initiate dictionary of dictionaries per cycle 
    dict_repeats = dict.fromkeys(cycles , {})
    dict_repeats['rgt'] = rgt #adds rgt to dictionary as useful metadata
    dict_repeats['rpt'] = rpt
    dict_repeats['lr'] = lr

    for cycle in cycles:
        
        file = glob.glob(os.path.join(atl06_dir, f'*processed_ATL06_*_{rgt}{cycle}*.h5'))[0]
        
        #Extract ATL06 data from the data dict per cycle, plus corresponding tide data 
        dict_repeats[cycle] = {'latitude': D[file]['latitude'],\
                               'h_li': D[file]['h_li'],\
                               'lat_corr': D[file]['lat_corr'],\
                               'lon_corr': D[file]['lon_corr'],\
                               'h_corr': D[file]['h_corr'],\
                               'atl06_quality_summary':D[file]['atl06_quality_summary'],\
                               'good_idx': D[file]['good_idx']}
        
        if tides==True:
            tidefile = glob.glob(os.path.join(atl06_dir, f'*ATL06_CATS2008_TIDES_*_{rgt}{cycle}*.h5'))[0]
            #Extract along-track tide and corresponding lat data 
            t = h5py.File(tidefile, 'r')
            tide = np.array(t['gt'+str(rpt)+lr]['land_ice_segments']['geophysical']['tide_ocean'])
            tide_lats = np.array(t['gt'+str(rpt)+lr]['land_ice_segments']['latitude'])

            #Set all values to nan where tide is set to max float value (i.e. nodata)
            if tide.max() > 50: 
                nanval = tide.max()
            else: 
                nanval = 9999999
            tide[tide == nanval] = np.nan
            tide_lats[tide == nanval] = np.nan

            tide_corr = tide[D[file]['good_idx']] #get tide values just at the latitudes where we have h_corr data
            tide_dict = {'tide_lats':tide_lats,\
                         'tide':tide,\
                         'tide_corr': tide_corr}
            dict_repeats[cycle].update(tide_dict)

    return dict_repeats


###############################################
## FUNCTIONS FOR EXTRACTING BED INFO ##########
###############################################

def find_nearest(array,value):
    
    '''Finds nearest point in an array to certain
    value (or points to values)'''
    
    idx = (abs(array-value)).argmin()

    return idx

def extract_bed_h(raster=None, A=None, B=None, np=1000, plot=False, ylabel='Bed Elevation [m]'):
    
    '''Extracts along-track bedrock elevation from input raster 
    (in this case the BedMachine bed geotiff, but can be used for
    any raster.)'''
    
    # Read raster data
    ds = gdal.Open(raster, 0)
    rb = ds.ReadAsArray()

    # Get start/end point of the map in lon/lat
    coords = ds.GetGeoTransform()

    # Get lon/lat and elevation
    lon = array([coords[0]+coords[1]*i for i in range(rb.shape[1])])
    lat = array([coords[3]+coords[5]*i for i in range(rb.shape[0])])
    elv = rb.flatten()

    # Make pairs of (lon,lat)
    points = array([array([x,y]) for y in lat for x in lon])

    # Define Profile using lon/lat
    ABx   = linspace(A[0],B[0],np)
    ABy   = linspace(A[1],B[1],np)
    profile = array([array([x,y]) for x,y in zip(ABx, ABy)])

    # Find Nearest points of profile in Map and Prepare for plot
    cross_h = []
    cross_lats = []
    
    #set nan value
    nanval = -9999

    for p in profile:
        lon_ind = find_nearest(lon,p[0])
        lat_ind = find_nearest(lat,p[1])
    
        x = lon[lon_ind]
        y = lat[lat_ind]
        h = rb[lat_ind][lon_ind] #Extract surface height at given lat/lon
        
        #only append non-nan values
        if h > nanval:
            cross_lats.append(y)
            cross_h.append(h)

    if plot == True:
        plt.figure()
        plt.xlabel('Latitude') 
        plt.ylabel(ylabel)
        plt.plot(cross_lats, cross_h)
    return cross_lats, cross_h

###############################################
## CALCULATING FITS OF SURFACE PROFILES #######
###############################################


def GenerateNewLats(dict_repeats, cycles, spacing=800): 
    
    '''Generate new lat spacing between min and max values
    across all cycles, with a spacing of 800 (currently)'''
    
    latmax = 0
    latmin = -180
    lenmin = 100000
    
    for cycle in cycles:
        if dict_repeats[cycle]['lat_corr'].max() < latmax:
            latmax = dict_repeats[cycle]['lat_corr'].max()
            idx = np.argmax(dict_repeats[cycle]['lat_corr'])
            lonmax = dict_repeats[cycle]['lon_corr'][idx] #get corresponding longitude at the max latitude (note - not necessarily the maximum lon)
        if dict_repeats[cycle]['lat_corr'].min() > latmin: 
            latmin = dict_repeats[cycle]['lat_corr'].min()
            idx = np.argmin(dict_repeats[cycle]['lat_corr'])
            lonmin = dict_repeats[cycle]['lon_corr'][idx] #get corresponding longitude at the min latitude (note - not necessarily the maximum lon)
        if dict_repeats[cycle]['lat_corr'].shape[0] < lenmin: 
            lenmin = dict_repeats[cycle]['lat_corr'].shape[0]

    if lenmin < spacing:
        spacing = lenmin 
    #print('Spacing: ', spacing)
    newlats = np.linspace(latmax, latmin, spacing)
    return newlats, lonmax, lonmin, latmax, latmin


def CalculateFits(dict_repeats, cycles, cycles_for_mean, s=0.7, tides=True, ref='mean',ref_cyc=None, extract_bed=False, bed=None):

    '''Function to calculate fits to each of the ATL06 elevation profiles per cycle. 
    Inputs: 
    - dict_repeats: output of ExtractData function 
    - cycles
    - cycles_for_mean: 
    - s: smoothing parameter (default = 0.7)
    - tides: set to True to also subtract along-track tides from ATL06 elevation and anomalies
    - ref: set to 'mean', 'neutral' or 'lowest'
    - ref_cyc: if using neutral or lowest tide, set to string of the reference tide e.g. '05' if that's the lowest tide 
    - extract_bed: set to True to extract along-track bed profile from raster
    - bed: file path to raster of bed elevation (must be in epsg 4326)'''
    
    #Create separate dictionary to be used for calculating mean profile
    fits = {}

    newlats, lonmax, lonmin, latmax, latmin = GenerateNewLats(dict_repeats, cycles, spacing=800) #change spacing as required
    dict_repeats['latmin'] = latmin
    dict_repeats['latmax']= latmax
    dict_repeats['lonmin'] = lonmin
    dict_repeats['lonmax']= lonmax

    # Iterate per cycle:
    for cycle in cycles:
        #Reverse array of lat and h_corr values because interpolation function cannot deal with decreasing x values (depends on whether track is ascending or descending)
        if dict_repeats[cycle]['lat_corr'][0]>dict_repeats[cycle]['lat_corr'][-1]:
            x = np.flip(dict_repeats[cycle]['lat_corr'])
            y = np.flip(dict_repeats[cycle]['h_corr'])
        else: 
            x = dict_repeats[cycle]['lat_corr']
            y = dict_repeats[cycle]['h_corr']

        #set weights (0 where nan, 1 if non-nan)
        w = np.where(np.isnan(y) == True, 0, 1)

        #1. find B-spline representation of the curve in a 2d plane. Outputs: 
            #tck = A 3-tuple (t, c, k) containing knot-points (t), coefficients (c) and the order (k) of the spline
            #fp = The weighted sum of squared residuals of the spline approximation (error metric)
            #ier = An integer flag about splrep success. Success is indicated if ier<=0.
            #msg = A message corresponding to the integer flag, ier.
        tck, fp, ier, msg= interpolate.splrep(x, y, w=w, s=s, full_output=1)

        #2. Evaluate the spline. Returns an array of values representing the spline function evaluated at the points in x.     
        y_ = interpolate.splev(newlats, tck, der=0) #can set ext=0 if want to extrapolate outside the given interval 
        dict_repeats[cycle]['y_'] = y_
        dict_repeats[cycle]['newlats_for_y_'] = newlats
        
        if cycle in cycles_for_mean:
            fits[cycle] = y_

        if ref != 'mean' and cycle == ref_cyc: #if using neutral / lowest tide, save these values for the reference track in dict_repeats
            ref_tck = tck
            ref_y_ = y_
            dict_repeats['ref_ATL06_y_'] = ref_y_
            dict_repeats['ref_ATL06_lats'] = newlats

        #3. Print out extra info  
        #print(f'Cycle {cycle}, s={s}: \n fp= {fp}\n ier= {ier}\n msg= {msg}')

    # a) Calculate mean profile between fits --> add to dictionary 
    mean_y_ = np.nanmean(np.array(list(fits.values())), axis=0)
    dict_repeats['mean_ATL06_y_'] = mean_y_
    dict_repeats['mean_ATL06_lats'] = newlats

    # b) Interpolate the mean y_ profile 
    mean_tck = interpolate.splrep(np.flip(dict_repeats[cycle]['newlats_for_y_']), np.flip(mean_y_), s=0)
    mean_profile_y_ = interpolate.splev(np.flip(newlats), mean_tck, der=0)

    # c) Calculate anomaly values between original h_corr values per cycle and the neutral tide fit
    for cycle in cycles:
        mean_y_ = interpolate.splev(dict_repeats[cycle]['lat_corr'], mean_tck, der=0)
        dict_repeats[cycle]['anom'] = dict_repeats[cycle]['h_corr'] - mean_y_
        if ref != 'mean':
            #Interpolate reference (neutr/lowest) tide heights at the individual cycle lat spacing
            ref_y_ = interpolate.splev(dict_repeats[cycle]['lat_corr'], ref_tck, der=0)
            dict_repeats[cycle]['anom_ref_tide'] = dict_repeats[cycle]['h_corr'] - ref_y_

    if tides == True:
        # d) Subtract tide amplitude from ATL06 elevations 
        for cycle in cycles: 
            dict_repeats[cycle]['hcorr_minus_tide'] = dict_repeats[cycle]['h_corr'] - dict_repeats[cycle]['tide_corr']    

        # e) Subtract tide amplitude from ATL06 elevation anomalies 
        for cycle in cycles: 
            dict_repeats[cycle]['anom_minus_tide'] = dict_repeats[cycle]['anom'] - dict_repeats[cycle]['tide_corr']
            if ref != 'mean': 
                dict_repeats[cycle]['anom_ref_minus_tide'] = dict_repeats[cycle]['anom_ref_tide'] - dict_repeats[cycle]['tide_corr']

    # f) Get bed profile and add to dictionary 
    if extract_bed == True:
        A = (lonmax, latmax)
        B=(lonmin, latmin)
        bed_lats, bed_h = extract_bed_h(bed, A=A, B=B, np=10000, plot=False)
        dict_repeats['bed']={'lats':bed_lats, 'h':bed_h}
    
    return dict_repeats

###############################################
## PLOTTING FUNCTIONS                       ###
## 1. Elevation surface profile             ###
## 2. Elevation anomaly                     ###
## 3. Combined figure (+ tides and bedrock) ###
###############################################

def PlotElevProfile(D, cycles, cycles_for_mean, track_dir, color_dict, xlim=None, ylim=None, plot_fits=True, plot_bed=True, plot_ref=True, save_fig=False, legend=True, figsize=[10,5],ref='mean',ref_cyc=None,extra_name=''): 
    
    '''Inputs: 
    D = dictionary (formerly dict_repeats)
    cycles = cycles to plot 
    cycles_for_mean = cycles used for fits to calculate the mean 
    plot_fits = set to True to plot the fits as well as the raw data
    save_fig = set to True to export the figure
    track_dir = path to directory for output per track
    figsize = set figure size [width, height]
    color_dict = dictionary to plot colors of each cycle
    
    Outputs: 
    '''

    bed_label = ''
    if plot_bed == True: 
        bed_label = '_bed'
        # Create two subplots and unpack the output array immediately
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=[10,6])

        #Plot elevation profile
        for cycle in cycles:
            #Plot fit and original data per cycle 
            if cycle in cycles_for_mean:
                ax1.plot(D['mean_ATL06_lats'], D[cycle]['y_'], label=f'Cycle {cycle} fit', color=color_dict[int(cycle)])
            latfilt = (D[cycle]['lat_corr'] > D['latmin']) & (D[cycle]['lat_corr'] < D['latmax']) #only plot within lat range covered by all tracks
            ax1.scatter(D[cycle]['lat_corr'][latfilt], D[cycle]['h_corr'][latfilt], label=f'Cycle {cycle} data', color=color_dict[int(cycle)], s=1)

        #plot reference elevation profile 
        if plot_ref == True:
            if ref =='mean': #If reference profile = mean of ATL06 fits:
                ax1.plot(D['mean_ATL06_lats'], D['mean_ATL06_y_'], label=f'Mean of ATL06 fits:\n{cycles_for_mean}', color='black')
            if ref =='neutral' or ref =='lowest': #If reference profile = cycle with most neutral tide
                ax1.plot(D['ref_ATL06_lats'], D['ref_ATL06_y_'], label=f'Mean of ATL06 fits:\n{cycles_for_mean}', color='black')

        #plot underlying bed profile (BedMachine)
        ax2.plot(D['bed']['lats'], D['bed']['h'], color='saddlebrown')
        ax2.fill_between(D['bed']['lats'], min(D['bed']['h']), D['bed']['h'], color='saddlebrown')

        #Formatting options
        plt.xlabel('Latitude')
        ax1.set_ylabel('ATL06 Surface Elevation (m)') 
        ax2.set_ylabel('BedMachine Bed Elevation (m)') 
        ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.suptitle(f'Track {D["rgt"]}, GT{D["rpt"]}{D["lr"].upper()}', weight='bold', size=16)
        plt.tight_layout();
        
    else: 
        plt.figure(figsize=figsize)
        for cycle in cycles:
            #Plot original data per cycle
            latfilt = (D[cycle]['lat_corr'] > D['latmin']) & (D[cycle]['lat_corr'] < D['latmax']) #only plot within lat range covered by all tracks
            plt.scatter(D[cycle]['lat_corr'][latfilt], D[cycle]['h_corr'][latfilt], label=f'Cycle {cycle} data', color=color_dict[int(cycle)], s=1)

            if plot_fits == True and cycle in cycles_for_mean:
                #Plot fit per cycle
                plt.plot(D[cycle]['newlats_for_y_'], D[cycle]['y_'], label=f'Cycle {cycle} fit', color=color_dict[int(cycle)])

        #Add the reference elevation profile (mean of ATL06 fits)
        if plot_ref == True:
            if ref =='mean': #If reference profile = mean of ATL06 fits:
                plt.plot(D['mean_ATL06_lats'], D['mean_ATL06_y_'], label=f'Mean of ATL06 fits:\n{cycles_for_mean}', color='black')
            if ref =='neutral' or ref =='lowest': #If reference profile = cycle with most neutral tide
                plt.plot(D['ref_ATL06_lats'], D['ref_ATL06_y_'], label=f'Mean of ATL06 fits:\n{cycles_for_mean}', color='black')
           
        #Formatting 
        plt.xlabel('Latitude')
        plt.ylabel('Elevation (m)')
        if legend== True:
            plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.title(f'Track {D["rgt"]}, GT{D["rpt"]}{D["lr"].upper()}', weight='bold', size=16)
        plt.tight_layout();
    if ylim is not None: 
        plt.ylim(ylim)
    if xlim is not None:
        plt.xlim(xlim)
    if save_fig == True:
        cyc_names = "".join(cycles)
        
        #Export figure
        outname_png = f'{track_dir}Track{D["rgt"]}_GT{D["rpt"]}{D["lr"].upper()}_elevation{bed_label}_cycs{cyc_names}_{ref}{extra_name}.png'
        outname_eps = f'{track_dir}Track{D["rgt"]}_GT{D["rpt"]}{D["lr"].upper()}_elevation{bed_label}_cycs{cyc_names}_{ref}{extra_name}.eps'
        plt.savefig(outname_png, dpi=300)
        plt.savefig(outname_eps, format = 'eps')

###############################################
###############################################
        
def PlotElevAnom(D, cycles, color_dict, track_dir, xlim=None, ylim=None,plot_bed=True, plot_tide=True, tide_position='all', save_fig=False, legend=True, plot_axhline=False, figsize=[10,5],ref='mean',ref_cyc=None, extra_name=''):

    tide_label = '_'
    bed_label = '_' 
    
    # Option to plot bedrock 
    if plot_bed == True:
        bed_label = '_bed'
        # Create two subplots and unpack the output array immediately
        f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=[10,8])
        for cycle in cycles:
            latfilt = (D[cycle]['lat_corr'] > D['latmin']) & (D[cycle]['lat_corr'] < D['latmax'])
            ax1.scatter(D[cycle]['lat_corr'][latfilt], D[cycle]['anom'][latfilt], label=f'Cycle {cycle}',color=color_dict[int(cycle)], s=2)
            if plot_tide == True:
                tide_label = '_tides'
                ax1.plot(D[cycle]['tide_lats'], D[cycle]['tide'], linestyle='dashed', label=f'Tide, Cycle {cycle}', color=color_dict[int(cycle)])
        
        #Plot bed profile in second subplot
        ax2.plot(D['bed']['lats'], D['bed']['h'], color='saddlebrown')
        ax2.fill_between(D['bed']['lats'], min(D['bed']['h']), D['bed']['h'], color='saddlebrown')
        
        #Formatting 
        plt.xlabel('Latitude')
        ax1.set_ylabel('Elevation Anomaly / Tide Amplitide (m)') 
        ax2.set_ylabel('BedMachine Bed Elevation (m)') 
        if legend== True:
            ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.suptitle(f'Track {D["rgt"]}, GT{D["rpt"]}{D["lr"].upper()}', weight='bold', size=16)
        plt.tight_layout()
        
    else: 
        plt.figure(figsize=figsize)
        for cycle in cycles:
            latfilt = (D[cycle]['lat_corr'] > D['latmin']) & (D[cycle]['lat_corr'] < D['latmax'])
            plt.scatter(D[cycle]['lat_corr'][latfilt], D[cycle]['anom'][latfilt], label=f'Cycle {cycle}', color=color_dict[int(cycle)], s=2)
            if plot_tide == True:
                tide_label = '_tides'
                #set filter to extract tide values for last 0.02 degrees lat along-track 
                filt =(D[cycle]['tide_lats']> (D['latmax']-0.02)) & (D[cycle]['tide_lats']< D['latmax'])
                
                #get latitudes for last 0.02 degrees along-track and reproject them further to the right (by 0.025 degrees lat) 
                tide_lats_plot = D[cycle]['tide_lats'][filt] + 0.025
                
                if ref == 'mean' or ref == 'neutral':
                    if tide_position == 'all':  #If want to plot all along-track tides 
                        
                        plt.plot(D[cycle]['tide_lats'], D[cycle]['tide'], linestyle='dashed', label=f'Tide, Cycle {cycle}', color=color_dict[int(cycle)])
                
                    if tide_position == 'overlap': #If only want to show tide at the edge of the figure (overlapping with anomaly data) 
                        plt.plot(D[cycle]['tide_lats'][filt], D[cycle]['tide'][filt], linestyle='dashed', label=f'Tide, Cycle {cycle}', color=color_dict[int(cycle)])
                
                    if tide_position == 'right': #If want to plot tide height to the right of the anomaly data 
                        plt.plot(tide_lats_plot, D[cycle]['tide'][filt], linestyle='dashed', label=f'Tide, Cycle {cycle}', color=color_dict[int(cycle)])
                
                if ref=='lowest': #when using the lowest tide as reference we convert tides to relative tides against lowest tide
                    try:
                        #1. calculate relative tide height against lowest tide
                        relative_tides = D[cycle]['tide']-D[ref_cyc]['tide'] 
                        
                        #2. Plot relative tide on the anomaly plots depending on position you want the dashed line (all, overlap or right)
                        if tide_position == 'all': #If want to plot all along-track tides 
                            plt.plot(D[cycle]['tide_lats'], relative_tides, linestyle='dashed', label=f'Tide, Cycle {cycle}', color=color_dict[int(cycle)])

                        if tide_position == 'overlap': #If only want to show tide at the edge of the figure (overlapping with anomaly data)
                            plt.plot(D[cycle]['tide_lats'][filt], relative_tides[filt], linestyle='dashed', label=f'Tide, Cycle {cycle}', color=color_dict[int(cycle)])

                        if tide_position == 'right': #If want to plot tide height to the right of the anomaly data 
                            plt.plot(tide_lats_plot, relative_tides[filt], linestyle='dashed', label=f'Tide, Cycle {cycle}', color=color_dict[int(cycle)])
              
                    except Exception as e: 
                        print('Cycle: ',cycle, e)
        
        # Plot line at y=0 to show the reference profile baseline 
        if plot_axhline==True: #If want to plot a reference line at y=0:
            if ref == 'mean':
                plt.axhline(y=0, xmax=1,label=f'Reference using cycle {ref_cyc} with {ref} tide', color='k')
            else:
                plt.axhline(y=0, xmax=1,label=f'Reference using cycle {ref_cyc} with {ref} tide', color=color_dict[int(ref_cyc)])
        
        plt.xlabel('Latitude')
        plt.ylabel('Elevation Anomaly / Tide Amplitide (m)') 
        if legend== True:
            plt.legend(title='Elevation anomaly:', bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.title(f'Track {D["rgt"]}, GT{D["rpt"]}{D["lr"].upper()}', weight='bold', size=16)
    
        plt.tight_layout()
    
    if ylim is not None: 
        plt.ylim(ylim)
    if xlim is not None:
        plt.xlim(xlim)
    
    #Export figure
    if save_fig == True:
        cyc_names = "".join(cycles)
        outname_png = f'{track_dir}Track{D["rgt"]}_GT{D["rpt"]}{D["lr"].upper()}_repeats_anom{tide_label}{bed_label}_cycs{cyc_names}_{ref}{extra_name}.png'

        outname_eps = f'{track_dir}Track{D["rgt"]}_GT{D["rpt"]}{D["lr"].upper()}_repeats_anom{tide_label}{bed_label}_cycs{cyc_names}_{ref}{extra_name}.eps'
        print(outname_png)
        plt.savefig(outname_png, dpi=300)
        plt.savefig(outname_eps, format = 'eps', dpi=300)
        
###############################################
###############################################
      