# agn_everest
Repository for testing [EVEREST 2.0](https://github.com/rodluger/everest) correction to K2 AGN and faint galaxies.

## Jupyter Notebooks
* everest_initial_test.ipynb - initial before/after comparison of EVEREST PLD 
* everest_manual_detrend.ipynb - applying EVEREST PLD to obejects that are not in EVEREST database
* C8_systematics_lc_templates.ipynb - Looking at the PSD of template light curves representing the K2 systematics we see
* C8_systematics_lc_templates_comparison.ipynb - calculating the X-correlation between observed/corrected light curves and their respective template systematic light curve

## python files
* analysis.py - contains lines importing useful modules and defines useful functions
* save_everest.py - script to manually detrend and save the light curves for objects in a given .csv file

## object_keys
Folder containing CSV files containing summary information of the objects we're testing
