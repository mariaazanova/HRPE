import numpy as np
import pandas as pd
import mne
import os
from collections import defaultdict
from support_functions import get_file_prefixes, setup_logging, process_files, compute_condition_differences, initialise_paths, clear_terminal, compute_condition_differences_update
import pickle 
import logging
import gc
import json

paths = initialise_paths()
root = paths['root']
evokedspath = paths['evokedspath']
tfrpath = paths['tfrpath']
GApath = paths['GApath']
pseudotrialpath = paths['pseudotrialpath']
logspath = paths['logspath']
derivativespath = paths['derivativespath']

logger = setup_logging(logspath, 'clearn14')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
perform_tfr_pseudotrials = True

subjects = get_file_prefixes(tfrpath,'_epo_tfr.hdf5','_epo_tfr.hdf5')

excel_filename = os.path.join(derivativespath, "num_ave_epochs.xlsx")
df_num_ave = pd.read_excel(excel_filename, index_col=0)

process_files(evokedspath, df_num_ave, pseudotrialpath, subjects, GApath)
compute_condition_differences_update(evokedspath, df_num_ave, subjects, GApath)

if perform_tfr_pseudotrials:
    process_files(tfrpath, df_num_ave, pseudotrialpath, subjects, GApath, is_tfr=True)
    compute_condition_differences_update(tfrpath, df_num_ave, subjects, GApath, is_tfr=True)
