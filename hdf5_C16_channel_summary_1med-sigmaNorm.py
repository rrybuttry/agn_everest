from analysis import *
import h5py
from astropy.convolution import convolve, Box1DKernel
import sys

if __name__=="__main__":

    if len(sys.argv)>1:
        channel = int(sys.argv[1])
    else:
        channel = 55

    campaign = 16
    mod,submod = get_submod(channel)
    submod-=1 # submods are 0-3 in the hdf5 file

    #hdf5_file = "/Users/rachelbuttry/K2/K2C%s_target_pixels.hdf5"%campaign
    hdf5_file = "/home/jackeline/Dropbox/Kepler/K2PixelMap_c16.hdf5"

    #hdf5_file = "/home/rachel/Research/K2/K2C%s_target_pixels.hdf5"%campaign
    # there are 3894 cadence points in c16
    cadences = 3894
    time = np.arange(cadences)/48.0

    with h5py.File(hdf5_file, 'r') as f:
        channel_epics = np.array(list(f['%s/%s/%s'%(campaign, mod,submod)].keys()))# look at objects in given channel).astype(int)

    # need the kepler magnitudes
    all_targets = pd.read_csv("/home/jackeline/Research/k2_pipeline/K2_official_targets/K2Campaign%stargets.csv"%campaign, skipinitialspace=True)

    # saving lcs to take the median
    arr = []

    with h5py.File(hdf5_file, 'r') as f:
        channel_epics = np.array(list(f['%s/%s/%s'%(campaign, mod, submod)].keys()))# look at objects in given channel).astype(int)

        rel_epics = channel_epics[np.isin(channel_epics, all_targets['EPIC ID'][np.logical_and(all_targets['magnitude'] > 13, all_targets['magnitude'] < 20)])]
        #print(len(rel_epics))

        # loop thru the objects
        for epic in rel_epics:
            d = np.array(f['%s/%s/%s/%s'%(campaign, mod, submod, epic)]['data'])
            lc_hdf5 = np.nansum(np.nansum(d, axis=1), axis=1)

            # handle spurious cadences
            lc = lk.LightCurve(time, flux=lc_hdf5)
            _, spurious_cad = lc.flatten().remove_outliers(return_mask=True) # remove spurious cadences
            lc_raw = lc.flux

            # interpolate
            cadno = np.arange(len(lc_raw)) # get an array to serve as our time/cadence measurement
            interped_vals = np.interp(cadno[spurious_cad], cadno[~spurious_cad], lc_raw[~spurious_cad])
            # replace spurious cadence values with the interpolated values
            lc_raw[spurious_cad] = interped_vals
            norm = np.std(lc_raw)
            lc_raw -= np.mean(lc_raw)
            lc_raw = lc_raw/norm

            smooth = convolve(lc_raw, Box1DKernel(350), boundary='extend')

            arr.append(lc_raw)

    #flux = np.nanmedian(np.atleast_2d(np.array(arr)), axis=0)
    #flux_smooth = convolve(flux, Box1DKernel(250), boundary='extend')

    #fig = plt.figure(figsize=(8,5))
    #plt.plot(time, flux[::-1], linewidth=4, alpha=0.8, color="#9ebcda")
    #plt.plot(time, flux_smooth[::-1], linewidth=10, color="#8856a7")
    
    #uncomment everything above to revert back 

    

    fig = plt.figure(figsize=(8,5))
    
    p = 50
    flux = np.percentile(np.atleast_2d(np.array(arr)),p, axis=0)
    flux_smooth = convolve(flux, Box1DKernel(250), boundary='extend')
    plt.plot(time, flux[::-1], linewidth=4, alpha=0.8, color="#9ebcda")
    plt.plot(time, (flux_smooth-flux_smooth.mean())[::-1], linewidth=10,  color="#8856a7")
    plt.ylim(-1.8,1.8)


    #plt.ylim(-.4,.4)
    plt.xticks([])
    plt.yticks([])
    plt.xlim(np.min(time), np.max(time))
    plt.tight_layout()
    plt.savefig("/home/jackeline/Research/k2_pipeline/submit_med/c16_channels/Channel%s.png"%channel, transparent=True, bbox_inches='tight', pad_inches=0)
