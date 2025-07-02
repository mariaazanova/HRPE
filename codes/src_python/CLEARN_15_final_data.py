import numpy as np
import pandas as pd
import mne
import os
import pickle 
from support_functions import get_file_prefixes, setup_logging, initialise_paths, find_best_match, extract_neural_data_multiple_roi, save_group_erp, save_channels_to_csv, clear_terminal, load_filename_mapping, list_mapped_files, handle_long_filename
import logging
import gc
import json

paths = initialise_paths()
root = paths['root']
epochspath = paths['epochspath']
evokedspath = paths['evokedspath']
tfrpath = paths['tfrpath']
GApath = paths['GApath']
pseudotrialpath = paths['pseudotrialpath']
derivativespath = paths['derivativespath']
logspath = paths['logspath']
avedfdifpath = paths['avedfdifpath']

logger = setup_logging(logspath, 'clearn15')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
baseline_window = (-0.2,0)
baseline_window_tfr = (-0.3,-0.2) # motivation: we're interested already in expectation-related erd
baseline_mode = 'zscore' # 'ratio','logratio','mean'
window_length = 0.05
# tfr settings
freqs = np.arange(0.1, 25, step=0.4)
n_cycles = freqs / 2.
decim = 4
freqs_of_interest = {} # change loading for freqs of interest from cfg
freqs_of_interest['high_alpha'] = [10, 13]
freqs_of_interest['low_alpha'] = [8, 10]
#freqs_of_interest['theta'] = [4, 7]
#freqs_of_interest['low_beta'] = [13, 16]
#freqs_of_interest['high_beta'] = [16, 25]
save_ecg = True
save_eog = True
save_breathing = False
eeg_of_interest = ['FCz', 'Pz']
ch_of_interest = ['FCz', 'Pz']
ch_ecg = ['EKG1']
ch_eog = ['EOGH','EOGV']
ch_resp = []
perform_tfr = True
save_GA_with_baseline_tfr = True 
save_individual_erps_with_baseline_tfr = True
save_GA_with_baseline = True 
save_individual_erps_with_baseline = False
# plotting options
save_channels_with_CI_based_on_ind_evo = True
save_channels_based_on_GA = True
only_topo_of_differences = False
only_tfr_differences = False

subjects = get_file_prefixes(epochspath,'.fif','.')

excel_filename = os.path.join(derivativespath, "num_ave_epochs.xlsx")
df_num_ave = pd.read_excel(excel_filename, index_col=0)

datapath = os.path.join(root, 'data.csv')
df = pd.read_csv(datapath)

avedatapath = os.path.join(root, 'average_data.csv')
avedf = pd.read_csv(avedatapath)

try:
    latencies = pd.read_csv(os.path.join(derivativespath, 'final_latencies.csv'))
except:
    latencies = pd.read_csv(os.path.join(derivativespath, 'all_significant_latencies.csv'))

rois = {}
if latencies.columns[0].startswith("Unnamed"):
    latencies = latencies.iloc[:, 1:]
for col in latencies.columns:
    val = latencies[col].iloc[0]  
    rois[col] = [round(float(val - window_length/2),3), round(float(val + window_length/2),3)]
rois['baseline'] = [-0.2,0] # baseline_window 
rois['tfr_baseline'] = [-0.3,-0.2] # baseline_window_tfr

if save_ecg:
    ch_of_interest = ch_of_interest + ch_ecg
if save_eog: 
    ch_of_interest = ch_of_interest + ch_eog
if save_breathing:
    ch_of_interest = ch_of_interest + ch_resp

# collect data from epochs to single-trial df
for subject in subjects:
    
    logger.info(f"Processing epochs of subject {subject}")

    pick_sub = (df.subID==int(subject))
        
    epochs = mne.read_epochs(os.path.join(epochspath, subject + '.fif'))  
    epochs.pick(ch_of_interest)

    extract_neural_data_multiple_roi(df, pick_sub, epochs, 
                                     ch_ecg, ch_eog, ch_resp, ch_of_interest, eeg_of_interest,
                                     save_ecg, save_eog, save_breathing, 
                                     rois, rename_roi="", 
                                     baseline_window=baseline_window, baseline_window_tfr=baseline_window_tfr)
    
    if save_GA_with_baseline or save_individual_erps_with_baseline:

        epochs.apply_baseline(baseline_window)

        extract_neural_data_multiple_roi(df, pick_sub, epochs, 
                                         ch_ecg, ch_eog, ch_resp, ch_of_interest, eeg_of_interest,
                                         save_ecg, save_eog, save_breathing, 
                                         rois, rename_roi="__withBase")
    del epochs

    if perform_tfr:

        logger.info(f"Processing tfr epochs of subject {subject}")
        
        epochs_tfr = mne.time_frequency.read_tfrs(os.path.join(tfrpath, subject +'_epo_tfr.hdf5'))
        epochs_tfr.pick(eeg_of_interest)

        extract_neural_data_multiple_roi(df, pick_sub, epochs_tfr, 
                                         ch_ecg, ch_eog, ch_resp, ch_of_interest, eeg_of_interest,
                                         save_ecg, save_eog, save_breathing, 
                                         rois, rename_roi="__TFR", 
                                         baseline_window=baseline_window, baseline_window_tfr=baseline_window_tfr,
                                         freqs_of_interest=freqs_of_interest)
        
        
        if save_GA_with_baseline_tfr or save_individual_erps_with_baseline_tfr:

            epochs_tfr.apply_baseline(baseline_window_tfr, mode=baseline_mode)
            
            extract_neural_data_multiple_roi(df, pick_sub, epochs_tfr, 
                                             ch_ecg, ch_eog, ch_resp, ch_of_interest, eeg_of_interest,
                                             save_ecg, save_eog, save_breathing, 
                                             rois, rename_roi="__TFR__withBase",
                                             baseline_window=baseline_window, baseline_window_tfr=None,
                                             freqs_of_interest=freqs_of_interest, dont_rescale=True)
        del epochs_tfr
    
    clear_terminal()
    gc.collect()

df.to_csv(datapath, index = False)
del df

# collect average subject-level values for difference between conditions 
dif_paths = [paths['evokedspath'], paths['tfrpath'], paths['pseudotrialpath']]
for dif_path in dif_paths:

    logger.info(f"Going through {dif_path}")

    mapping_file = os.path.join(dif_path, 'filename_mapping.csv')
    orig2short, short2orig = load_filename_mapping(mapping_file)
    files = list_mapped_files(dif_path)
    # Get lists of corrected and difference files
    virtual_files = [f['virtual'] for f in files]

    dif_files = [f for f in virtual_files if (f.startswith("dif_") or f.startswith("baseline_dif_") or f.startswith("long_dif_") 
                                                     or f.startswith("tfr_dif_") or f.startswith("tfr_baseline_dif_")
                                                    ) and f.endswith(".pkl")]

    for dif_file in dif_files:
        
        actual_file = orig2short.get(dif_file, dif_file)
        with open(os.path.join(dif_path, actual_file), 'rb') as f:
            dif_list = pickle.load(f)
        
        file_base = dif_file.replace(".pkl", "")

        if 'tfr' in dif_path:
            file_base = 'tfr_' + file_base
        if 'pseudo' in dif_path:
            file_base = 'pseudo_' + file_base
            if 'tfr_' in file_base:
                file_base = file_base.replace("tfr_", "")
                file_base = 'tfr_' + file_base

        if 'pseudo' not in dif_path:
            check_columns = find_best_match(file_base, list(df_num_ave)) 
            valid_subjects_mask = (df_num_ave[check_columns] > 0).all(axis=1).to_numpy()
            logger.info(f"For {file_base}, there are {sum(valid_subjects_mask)} valid subjects")
        else:
            valid_subjects_mask = np.full(len(df_num_ave), True, dtype=bool)
            logger.info(f"Processing {file_base}")

        temp_df = pd.DataFrame(index=avedf.index)

        subject_array = np.array(subjects)
        paired = [
            (sub, dif_el)
            for sub, dif_el, ok in zip(subject_array, dif_list, valid_subjects_mask)
            if ok
        ]
        logger.info(f"{len(paired)} valid (subject, dif) pairs after masking.")

        # Now iterate over the reduced list of (subject, dif) tuples
        for valid_subject_prefix, dif_list_element in paired:
            pick_sub = avedf['subID'] == int(valid_subject_prefix)

            rename_roi = "___" + file_base

            if 'tfr' in file_base:
                freqs_of_interest_current = freqs_of_interest
            else:
                freqs_of_interest_current = None
                
            if ('tfr' in file_base) and ('baseline' in file_base):
                dont_rescale_current = True
                baseline_window_tfr_current = None
                baseline_window_current = baseline_window
            elif 'baseline' in file_base:
                dont_rescale_current = False
                baseline_window_tfr_current = None
                baseline_window_current = None
            else:
                dont_rescale_current = False
                baseline_window_tfr_current = baseline_window_tfr
                baseline_window_current = baseline_window           

            extract_neural_data_multiple_roi(temp_df, pick_sub, dif_list_element, 
                                               ch_ecg, ch_eog, ch_resp, ch_of_interest, eeg_of_interest,
                                               save_ecg, save_eog, save_breathing, 
                                               rois, rename_roi, 
                                               baseline_window_current, baseline_window_tfr_current,
                                               freqs_of_interest_current, dont_rescale_current)
        
        tempdfpath = os.path.join(avedfdifpath, f"{file_base}.csv")
        tempdfpath = handle_long_filename(tempdfpath)
        temp_df.to_csv(tempdfpath, index = False)
        del temp_df

    clear_terminal()
    gc.collect()
del avedf  

# if we want, we can (again) collected info about GA and their CI from individual evokeds
if save_channels_with_CI_based_on_ind_evo:

    files = list_mapped_files(evokedspath)
    file_list = [f for f in files if f['virtual'].endswith(".pkl")]

    for file in file_list:

        file_actual = file['actual']
        file_virtual = file['virtual']

        full_path = os.path.join(evokedspath, file_actual)
        with open(full_path, 'rb') as f:
            data_list = pickle.load(f)  

        key = file_virtual.replace(".pkl", "")
        
        save_group_erp(key, data_list, GApath, save_ga=False, save_csv_with_ci=True, channels=ch_of_interest, return_ga=False, is_tfr=False)
    
clear_terminal()
gc.collect()

# we can also save just channels without CI from GA themselves, and also topograpies
files = list_mapped_files(GApath)
file_list = [f for f in files if f['virtual'].endswith(".fif")]
for file in file_list:

    file_actual = file['actual']
    file_virtual = file['virtual']

    full_path = os.path.join(GApath, file_actual)
    GA = mne.read_evokeds(full_path)[0]

    key = file_virtual.replace("_GA.fif", "")

    if save_channels_based_on_GA:

        save_channels_to_csv(key, GApath, GA, ch_of_interest)

    GA.pick('eeg')

    channel_names = GA.ch_names

    if only_topo_of_differences and 'dif_' not in key:
        continue

    if 'long_' in key:
        continue

    topo_df = pd.DataFrame({"Channel": channel_names})
    
    for roi, time_roi in rois.items():   
        
        avg_amplitudes = GA.get_data(tmin=time_roi[0], tmax=time_roi[1]).mean(axis=1)
        avg_amplitudes *= 1e6
        topo_df[roi] = avg_amplitudes

    topo_path = os.path.join(GApath, f"{key}_topo.csv")
    topo_path = handle_long_filename(topo_path)
    topo_df.to_csv(topo_path, index=False)
        
clear_terminal()
gc.collect()

# also save TFR at a channel and topographies
files = list_mapped_files(GApath)
file_list = [f for f in files if f['virtual'].endswith(".hdf5")]
for file in file_list:

    file_actual = file['actual']
    file_virtual = file['virtual']

    key = file_virtual.replace("_GA.hdf5", "")

    if only_topo_of_differences and only_tfr_differences and 'dif_' not in key:
        continue

    if not only_tfr_differences or (only_tfr_differences and 'dif_' in key):
        
        full_path = os.path.join(GApath, file_actual)
        GA = mne.time_frequency.read_tfrs(full_path)

        for channel in eeg_of_interest:
            
            tfr_df = pd.DataFrame(
                GA.get_data(picks=channel)[0],
                index=GA.freqs, 
                columns=GA.times
            )
            tfr_df.index.name = "Frequency (Hz)"
            tfr_df.columns.name = "Time (s)"
        
            tfr_file_name = f"{key}__{channel}.csv"
            tfr_path = os.path.join(GApath, tfr_file_name)
            tfr_path = handle_long_filename(tfr_path)
            tfr_df.to_csv(tfr_path)

    if not only_topo_of_differences or (only_topo_of_differences and 'dif_' in key):

        full_path = os.path.join(GApath, file_actual)
        GA = mne.time_frequency.read_tfrs(full_path)

        for freq, freq_band in freqs_of_interest.items():

            topo_df = pd.DataFrame({"Channel": channel_names})
        
            for roi, time_roi in rois.items(): 
                
                avg_amplitudes = GA.get_data(fmin=freq_band[0], fmax=freq_band[1], tmin=time_roi[0], tmax=time_roi[1]).mean(axis=1).mean(axis=1)
                if 'baseline' not in key:
                    avg_amplitudes *= 1e12                 
                topo_df[roi] = avg_amplitudes

            topo_path = os.path.join(GApath, f"{key}__{freq}_topo.csv")
            topo_path = handle_long_filename(topo_path)
            topo_df.to_csv(topo_path, index=False)  

    clear_terminal()
    gc.collect()
