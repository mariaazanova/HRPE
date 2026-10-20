import numpy as np
import pandas as pd
import mne
from mne import io 
import os
import re
from collections import defaultdict
from support_functions import get_file_prefixes, setup_logging, extract_evoked_latency, save_group_erp, initialise_paths, clear_terminal
from matplotlib import pyplot as plt
import gc
import logging
import json

paths = initialise_paths()
root = paths['root']
cleanpath = paths['cleanpath']
preprocpath = paths['preprocpath']
eventspath = paths['eventspath']
GApath = paths['GApath']
derivativespath = paths['derivativespath']
vispath = paths['cleaningvispath']
logspath = paths['logspath']

logger = setup_logging(logspath, 'clearn05')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
run_eeg_diagnostics = cfg["run_eeg_diagnostics"]
plot_each_sub = cfg["eeg_diag_plot_each_sub"]
save_ga = cfg["eeg_diag_save_ga"]
save_csv_with_ci = cfg["eeg_diag_save_csv_with_ci"]

subjects = get_file_prefixes(cleanpath,'.fdt','.')

if not run_eeg_diagnostics:
    logger.info('Not performing optional EEG diagnostics')
else:
    logger.info('Performing optional EEG diagnostics')

    files = os.listdir(preprocpath)
    pattern = re.compile(r'\d+_(\w+)\.set')
    keywords = set()
    for file in files:
        match = pattern.search(file)
        if match:
            keywords.add(match.group(1))
    if keywords:
        logger.info('Detected preprocessed files for diagnostics')
    
    ind_erp = defaultdict(list)
    group_erp = defaultdict(list)
    evoked_types = ['frn','p3a','p3b']
    latencies_of_interest = []
    
    for subject in subjects:
    
        subject_cleanpath= os.path.join(cleanpath, subject +'.set')
        subject_eventspath = os.path.join(eventspath, subject + '_allevents.xlsx')
    
        # load clean data after ica
        raw = io.read_raw_eeglab(subject_cleanpath, preload=True, eog=('EOGH','EOGV'))
        raw.set_channel_types({'EKG1':'ecg'})    
        raw.set_montage('standard_1005')    
    
        # load events from raw file and from excel 
        raw_events, _ = mne.events_from_annotations(raw, verbose=False)  
        events = pd.read_excel(subject_eventspath)
        # create a mapping (dict)
        event_dict = dict(zip(np.unique(events.type),np.unique(raw_events[:,2])))
        event_dict['Feedback/Correct/Approach'] = event_dict.pop('51')
        event_dict['Feedback/Correct/Avoid'] = event_dict.pop('52')
        event_dict['Feedback/Incorrect/Approach'] = event_dict.pop('53')
        event_dict['Feedback/Incorrect/Avoid'] = event_dict.pop('54')
        feedback_dict = {k: v for k, v in event_dict.items() if k.startswith('Feedback')}
    
        fig = raw.compute_psd(fmax=70).plot(show=False, sphere=(0.00, 0.02, 0.00, 0.09))
        fig.savefig(os.path.join(vispath, subject +'_clean_psd.jpg'), dpi=150, format='jpg')
        plt.close(fig)
    
        # epoch
        epochs = mne.Epochs(raw, raw_events, event_id=feedback_dict, tmin=-0.3, tmax=0.7, baseline=None,
                             event_repeated='drop', preload=True) 
        ind_erp['clean'].append(epochs.average(picks='all'))  
        del raw, epochs   
    
        for keyword in keywords: 
    
            subject_preprocpath= os.path.join(preprocpath, subject + '_' + keyword + '.set')
            
            # load preproc data before ica
            raw = io.read_raw_eeglab(subject_preprocpath, preload=True, eog=('EOGH','EOGV'))
            raw.set_channel_types({'EKG1':'ecg'})    
            raw.set_montage('standard_1005')    
        
            fig = raw.compute_psd(fmax=70).plot(show=False, sphere=(0.00, 0.02, 0.00, 0.09))
            fig.savefig(os.path.join(vispath, subject + '_' + keyword + '_psd.jpg'), dpi=150, format='jpg')
            plt.close(fig)
        
            # epoch
            epochs = mne.Epochs(raw, raw_events, event_id=feedback_dict, tmin=-0.3, tmax=0.7, baseline=None,
                                 event_repeated='drop', preload=True) 
            ind_erp[keyword].append(epochs.average(picks='all'))
            del raw, epochs   
    
        clear_terminal()

    for key, evoked_list in ind_erp.items():

        group_erp[key] = save_group_erp(key, evoked_list, GApath, save_ga=save_ga, save_csv_with_ci=save_csv_with_ci, channels=['FCz', 'Pz'], return_ga=True, is_tfr=False)
        
        if key=='clean':
            for evoked_type in evoked_types:
                latencies_of_interest.append(extract_evoked_latency(group_erp[key],evoked_type).item())
    
    gc.collect()

    for key in group_erp:
    
        fig = group_erp[key].plot(picks='FCz', show=False)
        fig.savefig(os.path.join(vispath, key +'_fcz.jpg'), dpi=150, format='jpg')
        plt.close(fig)
        fig = group_erp[key].plot(picks='Pz', show=False)
        fig.savefig(os.path.join(vispath, key +'_pz.jpg'), dpi=150, format='jpg')
        plt.close(fig)
        fig = group_erp[key].plot_topomap(times=latencies_of_interest, sphere=(0.00, 0.02, 0.00, 0.09), average = 0.05, show=False)
        fig.savefig(os.path.join(vispath, key +'_topo.jpg'), dpi=150, format='jpg')
        plt.close(fig)
        fig = group_erp[key].compute_psd(fmax=70).plot(show=False, sphere=(0.00, 0.02, 0.00, 0.09))
        fig.savefig(os.path.join(vispath, key +'_psd.jpg'), dpi=150, format='jpg')
        plt.close(fig)
    
        if plot_each_sub:
    
            for sub in range(len(subjects)):
    
                fig = ind_erp[key][sub].plot_topomap(times=latencies_of_interest, sphere=(0.00, 0.02, 0.00, 0.09), average = 0.05, show=False)
                fig.savefig(os.path.join(vispath, subjects[sub] + '_' + key +'_topo.jpg'), dpi=150, format='jpg')
                plt.close(fig)  
    
    del group_erp, ind_erp
    
    logger.info('Performed EEG diagnostics')
    clear_terminal()
