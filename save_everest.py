
# script to save everest de-trended shaya galaxies
# give it a list of EPIC IDs (input_csv has a column names "EPIC ID")

import everest
import pandas as pd
import numpy as np
import os

campaign = 8 
input_csv = "object_keys/c8_shaya_ch13_random50.csv"
folder = "../EVEREST_corrected"

def save_data(epic, campaign):
    """
    Function to save 4 numpy arrays into 2 files
    
    The first is a csv file with the following columns (NOT in time order, needs to be sorted):
        "time" : obj.time
        "flux" : obj.flux (not .fcor!)
        "flux_err" : (raw flux error ?)
    
    The second is a .npz containing
    "mask" : obj.mask (spurious cadence mask)
    
    """
    
    # download it if it's in EVEREST database, else detred maunally
    try:
        
        lc_everest = everest.Everest(epic, season=campaign)
        print("\t Found in EVEREST Database")
    except:
        print("\t Manually Detrending via nPLD")
        lc_everest = everest.detrender.nPLD(epic, season=campaign)
    

    # put it in a pandas dataframe
    df = pd.DataFrame(time = lc_everest.time,
        flux = lc_everest.flux,
        flux_err = lc_everest.frawerr)
    
    # save as a csv
    df.to_csv(folder+"%s_lc.csv"%epic, index=False)
    
    # save the spurious cadence mask
    np.save(folder+"%s_mask.npz"%epic, lc_everest.mask)
    
    

if __name__=="__main__":
    
    # create it if it doesn't exist
    if not os.path.exists(folder): os.mkdir(folder)
    
    # read in the input file
    df = pd.read_csv(input_csv)
    
    for epic in np.array(df["EPIC ID"]):
        print(epic)
        
        # check if files already exist
        if os.path.exists(folder+"%s_lc.csv"%epic) and os.path.exists(folder+"%s_mask.npz"%epic):
            print("\t Files already exist")
            pass
        else:
            print("\t Writing files")
            # run it if they don't exist
            save_data(epic, campaign)
            print("\t Done. ")