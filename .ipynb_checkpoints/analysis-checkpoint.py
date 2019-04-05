
# functions and fundamental values defined here!


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.stats import LombScargle

import lightkurve as lk
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
        t: time array
        y: flux array 
        
        f(optional): frequency range
    
    returns
        f: frequency range with first and last entry removed
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