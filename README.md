# agn_everest
[![DOI](https://zenodo.org/badge/178912855.svg)](https://zenodo.org/badge/latestdoi/178912855)
Repository for building off of [EVEREST 2.0](https://github.com/rodluger/everest) to create a correction to be applied to K2 AGN. Code is used in producing some figures in Moreno et al 2021

## Jupyter Notebooks
* arc_drift_obj.ipynb - Developing a way to determine/remove objects badly affected by arc drift
* everest_initial_test.ipynb - initial before/after comparison of EVEREST PLD 
* everest_manual_detrend.ipynb - applying EVEREST PLD to obejects that are not in EVEREST database
* everest_branch.ipynb - building off of EVEREST 2.0's PLD to create a correction method that won't remove variability
* C8_systematics_lc_templates.ipynb - Looking at the PSD of template light curves representing the K2 systematics we see
* C8_systematics_lc_templates_comparison.ipynb - calculating the X-correlation between observed/corrected light curves and their respective template systematic light curve
* PSD_slope_calc.ipynb - Creating and testing code to calculate slopes of PSDs in log-log space
* PSD_Combined.ipynb - Applying our corrections (see everest_branch.ipynb) and looking for relationships between fitted PSD slopes and physical parameters


## python files
* analysis.py - contains lines importing useful modules and defines useful functions
* save_everest.py - script to manually detrend and save the light curves for objects in a given .csv file

## object_keys
Folder containing CSV files containing summary information of the objects we're testing
