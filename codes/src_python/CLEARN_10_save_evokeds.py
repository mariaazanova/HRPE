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
evokedspath = paths['evokedspath']
GApath = paths['GApath']
derivativespath = paths['derivativespath']
logspath = paths['logspath']

logger = setup_logging(logspath, 'clearn10')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
save_GA_main = True # main conditions
save_GA_additional = True # additional control conditions
save_GA_sys_dia = False  # phase x miniblock conditions
# we need individual erps to plot averages after pseudotrial correction
save_individual_erps_main = True
save_individual_erps_additional = True
save_individual_erps_sys_dia = False
# in addition to saving evokeds without baseline (always), save evokeds with baseline applied before averaging
save_GA_with_baseline = True 
save_individual_erps_with_baseline = False
baseline_window = (-0.2,0)
# Save whole-epoch evokeds in addition to the short ones, they are never baseline corrected 
save_GA_long = False
save_individual_erps_long = False
get_difference = True # also save all pairwise differences of conditions
epoch_tmin_plot = -0.3
epoch_tmax_plot = 0.7
drop_bad_epochs = True

subjects = get_file_prefixes(epochspath,'.fif','.')

dfpath = os.path.join(root, 'data.csv') 
df = pd.read_csv(dfpath)

subtraction_order = {
    "outcome": [["Lose", "Win"]],          # Subtract "Lose" from "Win" for "outcome"
    "sysDiaFeedb": [["Diastole", "Systole"]] # Subtract "Diastole" from "Systole" for "sysDiaFeedb"
}
conditions, conditions_ga, conditions_ind = list_of_conditions(save_GA_main=save_GA_main, save_individual_erps_main=save_individual_erps_main, 
                                                               save_GA_additional=save_GA_additional, save_individual_erps_additional=save_individual_erps_additional,
                                                               save_GA_sys_dia=save_GA_sys_dia, save_individual_erps_sys_dia=save_individual_erps_sys_dia)

logger.info("Started collecting evokeds")

num_ave = defaultdict(list)

for condition in conditions:

    logger.info(f"Processing condition {condition}")

    ind_erp = defaultdict(list)
       
    for subject in subjects:

        logger.info(f"Processing condition {condition}, subject {subject}")
        
        epochs = mne.read_epochs(os.path.join(epochspath, subject + '.fif'))      
        
        if save_GA_long or save_individual_erps_long:
            if epochs.times[-1] == epoch_tmax_plot:
                logger.info("Epochs were already short. Not saving long evokeds.")
                save_GA_long = False
                save_individual_erps_long = False
            else:
                long_epochs = epochs.copy()
            
        if epochs.times[-1] > epoch_tmax_plot:
            logger.info("Initial epochs were long. Cropping.")
            epochs.crop(tmin = epoch_tmin_plot, tmax = epoch_tmax_plot)
        
        pick_sub = (df.subID==int(subject))
        sub_df = df.loc[pick_sub,].copy().reset_index(drop=True)   

        dict_evokeds, dict_num_ave = define_evokeds(epochs, sub_df, condition, drop_bad_epochs, False)
        if get_difference:
            dict_differences = calculate_evoked_differences(dict_evokeds, condition, subtraction_order)
            dict_evokeds.update(dict_differences) 
            del dict_differences
        for key in dict_evokeds:
            ind_erp[key].append(dict_evokeds[key])
            if not key.startswith('dif'):
                num_ave[key].append(dict_num_ave[key])
            
        if save_GA_with_baseline or save_individual_erps_with_baseline:
            dict_evokeds, dict_num_ave = define_evokeds(epochs, sub_df, condition, drop_bad_epochs, True, baseline_window) 
            if get_difference and save_individual_erps_with_baseline:
                dict_differences = calculate_evoked_differences(dict_evokeds, condition, subtraction_order)
                dict_evokeds.update(dict_differences)
                del dict_differences
            for key in dict_evokeds:
                new_key = f"{'baseline'}_{key}"
                ind_erp[new_key].append(dict_evokeds[key])

        del epochs

        if save_GA_long or save_individual_erps_long:
            dict_evokeds, dict_num_ave = define_evokeds(long_epochs, sub_df, condition, drop_bad_epochs, False)
            del long_epochs
            if get_difference and save_individual_erps_long:
                dict_differences = calculate_evoked_differences(dict_evokeds, condition, subtraction_order)
                dict_evokeds.update(dict_differences) 
                del dict_differences
            for key in dict_evokeds:
                new_key = f"{'long'}_{key}"
                ind_erp[new_key].append(dict_evokeds[key])   
        
        del dict_evokeds, dict_num_ave, sub_df

        logger.info(f"Processed condition {condition}, subject {subject}")
        gc.collect()
        clear_terminal()

    logger.info(f"Saving files for condition {condition}")
    for key in ind_erp.keys():

        save_this_key_ind = (key.startswith('baseline') and save_individual_erps_with_baseline) or (key.startswith('long') and save_individual_erps_long) or ((not key.startswith('long')) and (not key.startswith('baseline')))
        
        if save_this_key_ind and (condition in conditions_ind):
            filename = f"{key}.pkl"
            filepath = os.path.join(evokedspath, filename)
            with open(filepath, 'wb') as f:
                pickle.dump(ind_erp[key], f)

        save_this_key_GA = (key.startswith('baseline') and save_GA_with_baseline) or (key.startswith('long') and save_GA_long) or ((not key.startswith('long')) and (not key.startswith('baseline')))
            
        if save_this_key_GA and (condition in conditions_ga):

            save_group_erp(key, ind_erp[key], GApath, save_ga=True, save_csv_with_ci=True, channels=['FCz', 'Pz'], return_ga=False, is_tfr=False)

    logger.info(f"Finished condition {condition}")
    gc.collect()
    clear_terminal()
            
    del ind_erp
    
df_num_ave = pd.DataFrame(num_ave, index=subjects)
excel_filename = os.path.join(derivativespath, "num_ave_epochs.xlsx")
df_num_ave.to_excel(excel_filename)
del df_num_ave, num_ave