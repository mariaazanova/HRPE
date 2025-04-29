import numpy as np
import pandas as pd
import mne
from mne import io 
import os
from support_functions import get_file_prefixes, initialise_paths, setup_logging
import logging
import json

paths = initialise_paths()
root = paths['root']
cleanpath = paths['cleanpath']
preprocpath = paths['preprocpath']
eventspath = paths['eventspath']
GApath = paths['GApath']
derivativespath = paths['derivativespath']
vispath = paths['vispath']
logspath = paths['logspath']
epochspath = paths['epochs']

logger = setup_logging(logspath, 'clearn09')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
save_long_epochs = False
long_epoch_tmin = -3
long_epoch_tmax = 3
epoch_tmin_plot = -0.3
epoch_tmax_plot = 0.7
rejection_threshold = dict(eeg = 200e-6)
reject_tmin = 0.1
reject_tmax = 0.5
if save_long_epochs:
    epoch_tmin = long_epoch_tmin
    epoch_tmax = long_epoch_tmax
else:
    epoch_tmin = epoch_tmin_plot
    epoch_tmax = epoch_tmax_plot

subjects = get_file_prefixes(cleanpath,'.fdt','.')

dfpath = os.path.join(root, 'data.csv')  
df = pd.read_csv(dfpath)
df['bad_epoch'] = np.nan # we can combine it with bad_epochs identified in ecg, or keep separate columns

num_dropped_epochs = []

rpeaks = ['changeIBI4_3','changeIBI5_4','changeIBI6_5','changeIBI7_6', 'changeIBI3_2', 'changeIBI2_1', 'changeIBI8_7']

for rpeak in rpeaks:
    df[rpeak+'_2groups'] = None

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

    # epoch
    epochs = mne.Epochs(raw, raw_events, event_id=feedback_dict, tmin=epoch_tmin, tmax=epoch_tmax, baseline=None,
                         event_repeated='drop', preload=True) 
    epochs.save(os.path.join(epochspath, subject + '.fif'), overwrite=True, verbose=False)
    del raw

    pick_sub = (df.subID==int(subject))
    
    # bad epochs
    epochs_copy = epochs.copy()
    epochs_copy.reject_tmin = reject_tmin
    epochs_copy.reject_tmax = reject_tmax
    len_epochs_before = len(epochs)
    epochs_copy.drop_bad(reject=rejection_threshold)
    dropped_epochs = len_epochs_before - len(epochs_copy)
    num_dropped_epochs.append(dropped_epochs)
    df.loc[pick_sub, 'bad_epoch'] = (~np.in1d(epochs.selection, epochs_copy.selection))*1
    del epochs, epochs_copy

    for rpeak in rpeaks:
       df.loc[pick_sub, rpeak + '_2groups']= pd.qcut(df.loc[pick_sub, rpeak], 2, labels=['low','high']).values

# save df
df.to_csv(dfpath, index=False)
# save number of bad epochs
pd.DataFrame(num_dropped_epochs).to_excel(os.path.join(derivativespath, 'num_dropped_epochs.xlsx'))
# give a warning if someone has too many bad eeg epochs, dont drop for now
bad_epoch_thresh = 199
if np.any(np.array(num_dropped_epochs)>bad_epoch_thresh):
    print('Subject %s has too many bad epochs'%(np.array(subjects)[np.array(num_dropped_epochs)>bad_epoch_thresh]))
