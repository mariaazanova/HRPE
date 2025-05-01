import numpy as np
import pandas as pd
import mne
import os
from collections import defaultdict
from support_functions import get_file_prefixes, define_evokeds, calculate_evoked_differences, list_of_conditions, setup_logging, save_group_erp, initialise_paths, clear_terminal
import pickle
import gc
import logging
import json

paths = initialise_paths()
root = paths['root']
epochspath = paths['epochspath']
tfrpath = paths['tfrpath']
GApath = paths['GApath']
logspath = paths['logspath']

logger = setup_logging(logspath, 'clearn12')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
epoch_tmin_plot = -0.3
epoch_tmax_plot = 0.7
crop_longer_epochs = False
perform_tfr = True
subtract_overall_evoked = True
get_difference = True
drop_bad_epochs = True
hfreq = 0.1
lfreq = 25
step = 0.4
freqs = np.arange(hfreq, lfreq, step=step)
n_cycles = freqs / 2.
decim = 4
same_conditions_as_evo = False
if same_conditions_as_evo:
    print('load the same conditions') # load the same conditions for tfr as for evokeds
else: # compute only for main effects
    save_GA_main = True
    save_GA_additional = False
    save_GA_sys_dia = False
    save_individual_erps_main = True
    save_individual_erps_additional = False
    save_individual_erps_sys_dia = False
# baseline for erd ers 
save_GA_with_baseline_tfr = True 
save_individual_erps_with_baseline_tfr = True
baseline_window_tfr = (-0.3,-0.2) # motivation: we're interested already in expectation-related erd
baseline_mode = 'zscore'

subtraction_order = {
    "outcome": [["Lose", "Win"]],          # Subtract "Lose" from "Win" for "outcome"
    "sysDiaFeedb": [["Diastole", "Systole"]]
}

conditions, conditions_ga, conditions_ind = list_of_conditions(save_GA_main=save_GA_main, save_individual_erps_main=save_individual_erps_main, 
                                                               save_GA_additional=save_GA_additional, save_individual_erps_additional=save_individual_erps_additional,
                                                               save_GA_sys_dia=save_GA_sys_dia, save_individual_erps_sys_dia=save_individual_erps_sys_dia)

subjects = get_file_prefixes(epochspath,'.fif','.')

dfpath = os.path.join(root, 'data.csv') 
df = pd.read_csv(dfpath)

if not perform_tfr:
    logger.info("Skipping time frequency analysis")
else:
    logger.info("Started time frequency analysis")

    for i, subject in enumerate(subjects, start=1):
        logger.info(f"Started processing epochs for subject {i}/{len(subjects)}: {subject}")
        
        epochs = mne.read_epochs(os.path.join(epochspath, subject + '.fif'))    
        epochs.pick('eeg')
            
        # if epochs are too long initially, cut them from -0.5 to 0.9 and then again from -0.3 to 0.7 (to improve tfr at the edges)
        if epochs.times[-1] > epoch_tmax_plot:
            logger.info('Initial epochs were long. Cropping.')
            if crop_longer_epochs:
                epochs.crop(tmin = epoch_tmin_plot-0.2, tmax = epoch_tmax_plot+0.2)
                logger.info('Cropped epochs to a longer interval before tfr computation.')
            else:
                epochs.crop(tmin = epoch_tmin_plot, tmax = epoch_tmax_plot)
                logger.info('Cropped epochs to exact expected length.')
    
        if subtract_overall_evoked:
            epochs.subtract_evoked()   # clearly when we plot averages or do stats on averages, we need per-condition subtract_evoked. but that is not possible for single trial analysis on epochs + adds a lot of additional computation 
            
        epochs_tfr = epochs.compute_tfr(method='morlet', freqs=freqs, n_cycles=n_cycles, decim=decim, use_fft=True, n_jobs=-1, average=False, return_itc=False)        
        logger.info('Computed tfr.')
        del epochs
            
        if epochs_tfr.times[-1] > epoch_tmax_plot:
            epochs_tfr.crop(tmin = epoch_tmin_plot, tmax = epoch_tmax_plot)
            logger.info('Cropped epochs to exact expected length.')
            
        epochs_tfr.save(os.path.join(tfrpath, subject +'_epo_tfr.hdf5'), overwrite=True, verbose=False)
        del epochs_tfr
        logger.info(f"Finished processing epochs for subject {i}/{len(subjects)}: {subject}")
    logger.info("Finished processing epochs for all subjects")

    gc.collect()
    clear_terminal()

    logger.info("Started analysing conditions")
    for x, condition in enumerate(conditions, start=1):
        logger.info(f"Saving condition {x}/{len(conditions)}: {condition}")
        ind_erp = defaultdict(list)
           
        for i, subject in enumerate(subjects, start=1):
            logger.info(f"Processing subject {i}/{len(subjects)}: {subject}, condition {x}/{len(conditions)}: {condition}")
            
            epochs_tfr = mne.time_frequency.read_tfrs(os.path.join(tfrpath, subject +'_epo_tfr.hdf5'))      
            
            pick_sub = (df.subID==int(subject))
            sub_df = df.loc[pick_sub,].copy().reset_index(drop=True)   
    
            dict_evokeds, dict_num_ave = define_evokeds(epochs_tfr, sub_df, condition, drop_bad_epochs, False)
    
            if get_difference:
                dict_differences = calculate_evoked_differences(dict_evokeds, condition, subtraction_order)
                dict_evokeds.update(dict_differences) 
                del dict_differences
            for key in dict_evokeds:
                ind_erp[key].append(dict_evokeds[key])
    
            if save_GA_with_baseline_tfr or save_individual_erps_with_baseline_tfr:
                dict_evokeds, dict_num_ave = define_evokeds(epochs_tfr, sub_df, condition, drop_bad_epochs, True, baseline_window_tfr, baseline_mode) 
                if get_difference and save_individual_erps_with_baseline_tfr:
                    dict_differences = calculate_evoked_differences(dict_evokeds, condition, subtraction_order)
                    dict_evokeds.update(dict_differences)
                    del dict_differences
                for key in dict_evokeds:
                    new_key = f"{'baseline'}_{key}"
                    ind_erp[new_key].append(dict_evokeds[key])
    
            del dict_evokeds, dict_num_ave, sub_df, epochs_tfr  
            gc.collect()
            clear_terminal()

        logger.info(f"Finished collecting individual data for condition {x}/{len(conditions)}: {condition}")
    
        for key in ind_erp.keys():
    
            save_this_key_ind = (key.startswith('baseline') and save_individual_erps_with_baseline_tfr) or (not key.startswith('baseline'))
            
            if save_this_key_ind and (condition in conditions_ind):
                filename = f"{key}.pkl"
                filepath = os.path.join(tfrpath, filename)
                with open(filepath, 'wb') as f:
                    pickle.dump(ind_erp[key], f)
    
            save_this_key_GA = (key.startswith('baseline') and save_GA_with_baseline_tfr) or (not key.startswith('baseline'))
                
            if save_this_key_GA and (condition in conditions_ga):
                
                ind_erp[key] = [x for x in ind_erp[key] if x is not None]
                ga = mne.grand_average(ind_erp[key])
                ga.save(os.path.join(GApath, 'tfr_' + key + '_GA.hdf5'), overwrite=True, verbose=False)
                del ga

                save_group_erp(key, ind_erp[key], GApath, save_ga=True, return_ga=False, is_tfr=True)

        logger.info(f"Finished saving data for condition {x}/{len(conditions)}: {condition}")
                
        del ind_erp   
        gc.collect()
        clear_terminal()
