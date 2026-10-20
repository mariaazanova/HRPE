import numpy as np
import pandas as pd
import mne
from mne import io 
rng = np.random.default_rng(123)
import os
from collections import defaultdict
from support_functions import get_file_prefixes, setup_logging, get_miniblock_intervals, create_metadata_for_pseudotrials, get_ibi_group_thresholds, get_event_dict, define_evokeds, calculate_evoked_differences, load_evoked_response, list_of_conditions, initialise_paths, save_group_erp, clear_terminal
import pickle 
import logging
import gc
import json

paths = initialise_paths()
root = paths['root']
cleanpath = paths['cleanpath']
eventspath = paths['eventspath']
GApath = paths['GApath']
derivativespath = paths['derivativespath']
pseudotrialpath = paths['pseudotrialpath']
logspath = paths['logspath']

logger = setup_logging(logspath, 'clearn13')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
compute_all_pseudotrials = True 
# no matter which evokeds were computed, compute all possible pseudotrials, 
# if false, compute similar to saved ind erps (load them)
save_individual_erps_main = True
save_individual_erps_additional = False
save_individual_erps_sys_dia = False
max_retries = 500
num_sim = 10 # put 300 random feedback on the data 10 times 
split_IBI_into_groups_for_pseudotrials = 'all_ind_ibi_within_iteration' # 'slow_speed', 'all_ind_ibi', 'thresholds_from_data'
epoch_tmin_plot = -0.3
epoch_tmax_plot = 0.7
perform_tfr_pseudotrials = True # for all the same pseudotrials as defined in conditions, also do tfr analysis
subtract_overall_evoked = True
get_difference = True 
# tfr settings
hfreq = 0.1
lfreq = 25
step = 0.4
freqs = np.arange(hfreq, lfreq, step=step)
n_cycles = freqs / 2.
decim = 4

subtraction_order = {
    "outcome": [["Lose", "Win"]],          # Subtract "Lose" from "Win" for "outcome"
    "sysDiaFeedb": [["Diastole", "Systole"]]
}
conditions_pseudo = list_of_conditions(save_individual_erps_main=save_individual_erps_main,
                                       save_individual_erps_additional=save_individual_erps_additional,
                                       save_individual_erps_sys_dia=save_individual_erps_sys_dia,
                                       pseudotrial_conditions=True, compute_all_pseudotrials=compute_all_pseudotrials)

subjects = get_file_prefixes(cleanpath,'.fdt','.')

logger.info("Started simulations")
for subject in subjects:

    logger.info(f"Processing subject {subject}")

    subject_cleanpath= os.path.join(cleanpath, subject +'.set')
    subject_eventspath = os.path.join(eventspath, subject + '_allevents.xlsx')
     
    # load events, compute IBIs
    events = pd.read_excel(subject_eventspath)
    events.drop('duration', axis=1, inplace=True)
    feedb_latencies = events.loc[((events.type=='51')|(events.type=='52')|
                                  (events.type=='53')|(events.type=='54')),'latency'].values
    time_range_blocks = [events.loc[events.type=='101','latency'].values, 
                         events.loc[events.type=='102','latency'].values, 
                         events.loc[events.type=='103','latency'].values]
    Rpeak_latencies = events.loc[(events.type=='Rpeak'),'latency'].values
    IBIs = np.diff(Rpeak_latencies)
    change_IBI = np.diff(IBIs)
    IBIs = IBIs[1:]
    Rpeak_latencies = Rpeak_latencies[2:]
    
    # load raw data
    raw = io.read_raw_eeglab(subject_cleanpath, preload=True, eog=('EOGH','EOGV'))
    raw.set_channel_types({'EKG1':'ecg'})    
    raw.set_montage('standard_1005')    
    sfreq = int(raw.info['sfreq'])

    miniblock_intervals = get_miniblock_intervals(feedb_latencies, time_range_blocks, sfreq)
    thresh_ibi, thresh_change_ibi = get_ibi_group_thresholds(split_IBI_into_groups_for_pseudotrials)

    sub_sim = defaultdict(list) # all simulated ERPs of this subject

    for i in range(num_sim):

        retries = 0
        successful_simulation = False
        
        while not successful_simulation and retries < max_retries:
            retries += 1
            logger.info(f"  Simulation {i + 1}/{num_sim}, Attempt {retries}/{max_retries}")

            new_feedback_latencies = []
            for time_range in time_range_blocks:
                b = time_range[1]-sfreq
                a = time_range[0]+sfreq
                new_feedback_latencies.extend(sorted(rng.integers(low=a, high=b, size=100)))   
        
            meta = create_metadata_for_pseudotrials(new_feedback_latencies, Rpeak_latencies, IBIs, change_IBI, 
                                                    miniblock_intervals, thresh_ibi, thresh_change_ibi)
            new_events, new_event_dict = get_event_dict(new_feedback_latencies)
    
            epochs = mne.Epochs(raw, new_events, event_id=new_event_dict, tmin=epoch_tmin_plot, tmax=epoch_tmax_plot, baseline=None, event_repeated='drop', preload=True)
            logger.info(f"    Created {len(epochs)} epochs (expected: 300).")

            if len(epochs) < 300:
                logger.warning(f"Simulation {i + 1}, Attempt {retries}: Created {len(epochs)} epochs. Retrying...")
                continue  
            if len(meta) != len(epochs.events):
                logger.warning(f"Simulation {i + 1}, Attempt {retries}: Metadata length {len(meta)} "
                                f"does not match epochs {len(epochs.events)}. Retrying...")
                continue

            successful_simulation = True

        if not successful_simulation:
            logger.error(f"Simulation {i + 1} failed after {max_retries} retries. Skipping simulation.")
            continue
        
        if perform_tfr_pseudotrials:
            if subtract_overall_evoked:
                tfr_epochs = epochs.copy().subtract_evoked() # subtract evoked could be done per condition (more precise) - comparable to code 12, tfr
            else: 
                 tfr_epochs = epochs.copy()           
            tfr_epochs = tfr_epochs.compute_tfr(method='morlet', freqs=freqs, n_cycles=n_cycles, decim=decim, use_fft=True, n_jobs=-1, average=False, return_itc=False)        

        for condition in conditions_pseudo:
            # if we want baseline-corrected pseudotrials, add baseline within define_evokeds
            dict_evokeds, dict_num_ave = define_evokeds(epochs, meta, condition)
            for key, evoked in dict_evokeds.items():
                sub_sim[key].append(evoked)
                logger.info(f"  Key: {key}, epochs in simulation: {dict_num_ave[key]}")
            del dict_evokeds, dict_num_ave
       
            if perform_tfr_pseudotrials:
                # if we want baseline-corrected pseudotrials, add baseline within define_evokeds
                dict_evokeds, dict_num_ave = define_evokeds(tfr_epochs, meta, condition)
                for key, evoked in dict_evokeds.items():
                    sub_sim['tfr_'+key].append(evoked)
                del dict_evokeds, dict_num_ave

        del epochs

    del raw

    for key, evoked_list in sub_sim.items():
        
        evoked_list = [x for x in evoked_list if x is not None]
        count = len(evoked_list)
        logger.info(f"  Key: {key}, Successful simulations before grand average: {count}")
        ga = mne.grand_average(evoked_list)        
        if key.startswith('tfr'):
            ga.save(os.path.join(pseudotrialpath, subject + '_' + key +'_ave.hdf5'), overwrite=True, verbose=False)
        else:
            ga.save(os.path.join(pseudotrialpath, subject + '_' + key +'_ave.fif'), overwrite=True, verbose=False)
        del ga
    
    del sub_sim
    
    clear_terminal()
    gc.collect()

    files = [f for f in os.listdir(pseudotrialpath) if f.endswith("_ave.hdf5") or f.endswith("_ave.fif")]
unique_vars = list(set(string for sublist in conditions_pseudo for string in sublist))

keys = set()
for f in files:
    if any(f.endswith(ext) for ext in ["_ave.hdf5", "_ave.fif"]):
        parts = f.split('_', 1)  
        if len(parts) == 2:  
            remaining = parts[1]  
            key = remaining.rsplit('_', 1)[0]  
            keys.add(key)

if get_difference:
    logger.info("Started collecting differences")
    
    for i, subject in enumerate(subjects, start=1):
        logger.info(f"Processing subject {i}/{len(subjects)}: {subject}")
        dict_evokeds = {}
    
        for key in keys:
            evoked = load_evoked_response(subject, key, files, pseudotrialpath)
            dict_evokeds[key] = evoked
    
        dict_differences = calculate_evoked_differences(dict_evokeds, unique_vars, subtraction_order)
        del dict_evokeds
        if not dict_differences:  
            logger.warning(f"No differences computed for subject {subject}. Skipping subject.")
            continue  

        for key, evo in dict_differences.items():

            difdictpath = os.path.join(pseudotrialpath, f"{key}.pkl")
            
            if i>1:
                with open(difdictpath, 'rb') as f:
                    dif_list = pickle.load(f)
                dif_list.append(evo)
                with open(difdictpath, "wb") as f:
                    pickle.dump(dif_list, f)

                if i==len(subjects):
                    try:
                        if key.startswith('tfr_'): 
                            key_for_save_group_erp = "pseudo_" + key[4:]
                            save_group_erp(key_for_save_group_erp, dif_list, GApath, save_ga=True, return_ga=False, is_tfr=True)
                        else:
                            key_for_save_group_erp = "pseudo_" + key
                            save_group_erp(key_for_save_group_erp, dif_list, GApath, save_ga=True, save_csv_with_ci=True, channels=['FCz', 'Pz'], return_ga=False, is_tfr=False)            
                    except Exception as e:
                        logger.error(f"Error saving grand average for {key}: {e}")
                    
                del dif_list, evo
            
            else:
                dif_list = [evo] 
                with open(difdictpath, "wb") as f:
                    pickle.dump(dif_list, f)
                del dif_list, evo

        del dict_differences
        gc.collect()
        clear_terminal()

logger.info("Started GA computation loop")

for key in keys:
    logger.info(f"Processing key: {key}")

    evoked_list = []

    for subject in subjects:
        evoked = load_evoked_response(subject, key, files, pseudotrialpath)
        evoked_list.append(evoked)   
        del evoked
        
    try:
        if key.startswith('tfr_'): 
            key_for_save_group_erp = "pseudo_" + key[4:]
            save_group_erp(key_for_save_group_erp, evoked_list, GApath, save_ga=True, return_ga=False, is_tfr=True)
        else:
            key_for_save_group_erp = "pseudo_" + key
            save_group_erp(key_for_save_group_erp, evoked_list, GApath, save_ga=True, save_csv_with_ci=True, channels=['FCz', 'Pz'], return_ga=False, is_tfr=False)            
    except Exception as e:
        logger.error(f"Error saving grand average for {key}: {e}")

    del evoked_list
    gc.collect()
    clear_terminal()