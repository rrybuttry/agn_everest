
# functions and fundamental values defined here!
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.stats import LombScargle
import lightkurve as lk

from astropy.convolution import convolve, Box1DKernel # for boxcar smoothing 

from lightkurve.correctors import PLDCorrector
import everest 

import logging
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)


#-------------------------------------------Extracting correct K2 light curve fits file------------------------------

def mast_path(epic, c, dtype='lightcurves'):
    """
    Returns file directory and filename of K2 light curve fits files on MAST server
    (needs to be edited)
    """
    if len(str(c)) < 2: c_str = '0'+str(c)
    else: c_str = str(c)
    
    epic = str(int(epic))
    XXXX = epic[:4]
    YY = epic[4:6]
    #ZZZ = epic[6:]
    
    
    #dtype = 
    dir_path = 'https://archive.stsci.edu/pub/k2/lightcurves/c'+c+'/'+XXXX+'00000/'+YY+'000'
    fits_file = 'ktwo%s-c%s_llc.fits'%(epic, c_str)
    return dir_path, fname
    
def local_path(epic, c):
    """
    Returns file directory and filename of K2 target pixel files
    """
    if len(str(c)) < 2: c_str = '0'+str(c)
    else: c_str = str(c)
    # get path to fits
    #mastDownload_path = '/home/rachel/Research/K2/mastDownload/'
    #fpath = mastDownload_path+'K2/ktwo%s-c%s_lc/ktwo%s-c%s_lpd-targ.fits.gz'%(epic, c_str, epic, c_str)
    #dir_path = mastDownload_path+'K2/ktwo%s-c%s_lc/'%(epic, c_str)

    #dir_path = '/home/rachel/.lightkurve-cache/mastDownload/K2/ktwo%s-c%s_lc/'%(epic, c_str)
    dir_path = '~/.lightkurve-cache/mastDownload/K2/ktwo%s-c%s_lc/'%(epic, c_str)
    fname = 'ktwo%s-c%s_lpd-targ.fits.gz'%(epic, c_str)
    return dir_path, fname


def get_submod(channel):
    """
    function to convert K2 channel number to module
    
    returns
        mod: 2-24
        submod: 1-4
    """

    #if channel > 85 or channel < 1:
    #    raise ValueError('Invalid Channel Number: %i'%channel)

    offset = 2
    if channel%4==0:
        offset -= 1
    if channel > 12:
        offset += 1
    if channel > 72:
        offset += 1

    mod = channel//4 + offset
    
    submod = channel-(mod-offset)*4
    if submod==0: submod=4
    
    return mod, submod

#-------------------------------------------Plotting power spectrum density------------------------------------------

# K2 frequncy range we're probing 
#   80 day campaigns, minumum 30 minute spacings
#   half-hours = days*24*2
#   sec = min*60
k2_freq = np.fft.rfftfreq(80*48, 30.0*60)


def LS_PSD(t,y, f=None):
    """
    Normalized Lomb Scargle Power Spectrum Density
    
    args
        t : time array
        y : flux array 
        
        f(optional) : frequency range
    
    returns
        f : frequency range with first and last entry removed
        power_ls: lomb-scargle power
    
    """
    N = len(t) # number of datapoints    
    #f = np.fft.rfftfreq(N**2, (t[-1] - t[0])) # frequency range (from Uttley et al. 2002)
    if f is None:
        f = np.fft.rfftfreq(N, (t[1] - t[0])) # frequency range

    # Compute the LS based power spectrum estimates
    model = LombScargle(t, y)
    power_ls = model.power(f[1:-1], method="fast", normalization="psd")

    # >>> To get the LS based PSD in the correct units, normalize by N <<<
    power_ls /= N
    
    return f[1:-1], power_ls

def plot_lc_PSD(time, flux, ax1,ax2,l="", f=None, **kwargs):
    """
    Plots lightcurve and PSD in ax1 and ax2 respectively
    
    args
        time: time array
        flux: flux array
        ax1: matplotlib axis (subplot) for flux v. time
        ax2: matplotlib axis (subplot) for power v. freq
        l: (optional) label for input data (e.g. raw, corrected)
    
    """
    
    # sec = days*86400
    f,power = LS_PSD(time*86400,flux, f=f)
    #f,power = powerSpectrum(time,flux)
    
    
    ax1.plot(time,flux, label=l+" light curve", **kwargs)
    ax1.set_xlabel("Time - 2454833[BKJD days]")
    ax1.set_ylabel("Flux")
    ax1.legend()
    
    ax2.plot(f,power, label=l+" PSD", **kwargs)
    ax2.set_xlabel("frequency [Hz]")
    ax2.set_ylabel("power [$\mathrm{ppm}^2/\mathrm{Hz}$]")
    ax2.set_yscale("log")
    ax2.set_xscale("log")  
    ax2.legend()
    
#-------------------------------------------Statistics and metrics------------------------------------------

def norm_xcorr(a, v):
    """
    Compute the normalized cross-correlation (lag plot)
    
    args
        a : time series #1
        v : time series #2
        
    returns
        normlized cross correlation (where 1.0 is 100% correlated)
    """
    a = (a - np.mean(a)) / (np.std(a) * len(a))
    v = (v - np.mean(v)) /  np.std(v)

    return np.correlate(a, v, mode='full') # extracta the lag times


def gini(array):
    """
    Function to calaculate the Gini coefficient
    (code taken from: https://github.com/oliviaguest/gini)
    """
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    
    # All values are treated equally, arrays must be 1d:
    array = array.flatten()

    # remove NaNs
    array = array[~np.isnan(array)]

    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array += 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1,array.shape[0]+1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array)))


#-------------------------------------------Data pre-processing------------------------------------------

# (boolean) mask to apply to a K2 capaign 8 light curve to make it the same shape as the lightkurve time series
#  3852 cadence points in a complete 80 day lc
_test_lc_c8 = lk.search_targetpixelfile(220292052,campaign=8).download().to_lightcurve(aperture_mask='all')
c8_lk_mask = np.isin(np.arange(3852+1), _test_lc_c8.cadenceno-_test_lc_c8.cadenceno[0])

def mask_spurious(time,flux):
    """
    Function to apply lightkurve.flatten method of calculating spurious cadences
    
    returns
        spurious_cadences: boolean mask 
    """
    lc = lk.LightCurve(time=time, flux=flux)
    _, spurious_cadences = lc.flatten().remove_outliers(return_mask=True)
    
    return spurious_cadences

def ind_to_bool(range_arr, mask_in):
    """
    EVEREST K2 masks are arrays containing the  indicies of the relevant cadences.
    This function converts t
    (may remove this function)
    
    args
        range_arr: an np.arange object that represents the ordered 
            indices of the array the mask will be applied to
        mask_in: np array containing the indices of the mask(s)
        
    returns
        a boolean mask 
        
    
    """
    
    return np.isin(range_arr, mask_in)

def sort_by_time(t,f,mask=None):
    """
    Fucntion to sort time series by time stamp/
    """

    sorting = np.array(t).argsort()
    t_sorted = t[sorting]
    f_sorted = f[sorting]
        
    return t_sorted, f_sorted

def interp_missing(c, f, gap=None):
    """
    Function to "interpolate" all the missing cadence points
    
    args
        c : cadence number (must be in ints(?))
        f : flux at given cadence 
        
    returns
        c : complete cadence numbers
        f : flux at the complete cadence numbers
    """

    missing_c = np.setdiff1d(range(min(c), max(c)+1), c)

    # incorporate mid campaign gap if necessary
    if gap is not None: missing_c = missing_c[np.logical_not(np.isin(missing_c, gap))]


    for cm in missing_c:
        # get index right after "missing" time value
        ind = np.argwhere(c > cm)[0]
        # interpolate the "missing" corrected flux values (take the average--linear interp)
        missing_f = (f[ind-1] + f[ind])/2.0

        # insert them into the correct locations in the arrays
        f = np.insert(f, ind, missing_f)
        c = np.insert(c, ind, cm)

    return c, f

def bin_smoothing(x,y, N=30):
    """
    Function to perform smoothing on a dataset 
    Computes the average of bins for data binned with N data points.
    (Note: takes both x and y to make my life easier)
    
    args
        x : array of x dimension of the data
        y : array of y dimension of the data
        N : number of data points in bin (default N = 30)
    returns
        x_smooth : smoothed x-array
        y_smooth : smoothed y-array
    """
    x_smooth = []
    y_smooth = []

    i = 0
    while i < len(x):
        if i+N > len(x):
            x_smooth.append(np.average(x[i:]))
            y_smooth.append(np.average(y[i:]))
        else:
            x_smooth.append(np.average(x[i:i+N]))
            y_smooth.append(np.average(y[i:i+N]))
        i += N
        
    return np.array(x_smooth), np.array(y_smooth)

def smooth_fit(x, y, smoothing='boxcar', N=30):
    """
    Function to smooth and perform linear fit to x,y data
    
    args
        x : array of x dimension of the data
        y : array of y dimension of the data
        smoothing : method to smooth the data
            'boxcar' : (default) smoothed data is same size as input data
            'bin' : (see bin_smoothing) smoothed data is smaller shape than input data

    returns
        m : slope of smoothed fit
        b : y-intercept of smoothed fit
    """
    
    # smooth
    if smoothing=='boxcar':
        # boxcar smoothing
        y_smooth = convolve(y, Box1DKernel(N), boundary='extend')
        x_smooth = x
    elif smoothing=='bin':
        bin_smoothing(x,y,N)
    else:
        raise ValueError("Invalid smoothing type specified '%s'"%smoothing)
        
        
    # fit 
    m, b = np.polyfit(x_smooth, y_smooth, 1)
    
    return m,b



#-------------------------------------------K2 EVEREST PLD based correction------------------------------------------

def K2_correction(epic, campaign=None, lttr='cbv', runpld=True):
    """
    performs a PLD based correction on K2 objects
    
    args
        epic : EPIC ID (int)
        campaign : K2 campaign (int)
        runpld(optional) : if true, nPLD will be manually run when the EPIC ID is not found in the EVEREST database
        lttr(optional) : long-timescale trend removal method 
            'cbv' subtracts 2 cotrending basis vectors
            'linear' subtracts linear slope
        
    returns
        time : array of timestamps for observations (in days)
        flux_corrected : array of corrected flux for object
    """  

    try:
        # get an everest light curve 

        # -----------------get EVEREST PLD-----------------
        lc = everest.Everest(epic, season=campaign, mission='k2')

    except:
        # there is no corrected light curve in the EVEREST database
        logging.info("No EVEREST light curve found")
        
        if runpld:
            logging.info("Running EVEREST nPLD (this may take a while)")
        
            # -----------------run EVEREST nPLD-----------------
            lc = everest.nPLD(epic, season=campaign, mission='k2')
            logging.info("Done")
            
        else:
            logging.info("Exiting funtion returning None; run again and set 'runpld=True' to run the EVEREST correction")
            return None
    
    if lttr=='cbv':
        # calculate correction based on (2) cotrending basis vectors
        lc.cbv_num = 2
        lc.compute()
        
        
        # put flux/cadences into an array
        # (there are 3852 cadences in a given 80 day campaign)
        cad = np.arange(3852+1)
        flux_pld = lc.flux

        # turning indices found to be "bad" into a boolen mask to apply
        mask = (np.isin(cad, np.concatenate([lc.nanmask, lc.badmask, lc.mask])))

        # interpolate the spurious cadences
        interped_vals = np.interp(cad[mask], cad[~mask], flux_pld[~mask])
        # replace spurious cadence values with the interpolated values
        flux_pld[mask] = interped_vals

        
        return lc.time, flux_pld
        
    elif lttr=='linear':

        # -----------------do addtional cut and slope subtraction-----------------
        # put flux/cadences into an array
        # (there are 3852 cadences in a given 80 day campaign)
        cad = np.arange(3852+1)
        flux_pld = lc.flux

        # turning indices found to be "bad" into a boolen mask to apply
        mask = (np.isin(cad, np.concatenate([lc.nanmask, lc.badmask, lc.mask])))

        # interpolate the spurious cadences
        interped_vals = np.interp(cad[mask], cad[~mask], flux_pld[~mask])
        # replace spurious cadence values with the interpolated values
        flux_pld[mask] = interped_vals

        # 30 mintue intervals between cadences 
        cutoff_day = 3*24*2
        #cutoff = np.logical_and(cad>cutoff_day, cad<cad[-1]-5*cutoff_day)
        cutoff = cad>cutoff_day
        # finding linear fit 
        m,b = np.polyfit(cad[cutoff], flux_pld[cutoff], 1)

        # subtracting it
        flux_corrected = flux_pld[cutoff] - (m*cad[cutoff])

        # time from the light curve
        time = lc.time[cutoff]  

        return time, flux_corrected
    
    else:
        raise ValueError("Invalid method in lttr arg. Use either 'cbv' or 'linear'")


    