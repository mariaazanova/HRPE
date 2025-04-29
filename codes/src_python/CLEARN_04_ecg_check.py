# this is a placeholder for a code that:
# (1) performs identification of tpeaks using trapeziums approach, if required
# (2) plots subject-level trials and corresponding R-peak, T-peak events for manual inspection and adjustment

# for each subject, plot 15 figures (each with 20 trials) 
#  -5 before to 5 sec after, plot ecg and indentified peaks
# raw samples are shown on x axis of figures straight away
# we can mark a trial as bad or replace latency
# then a file is saved that contains all epochs of all subs and info if they are good based on ecg
# we also replace event latencies for rpeaks and t peaks in the original events tables
# we can also remove R-peaks or T-peaks from events or add them, not only correct
# automatic quality check: check that there are at least 150 good ecg epochs per sub, only give warn, dont remove at this stage

import numpy as np
import pandas as pd
import mne
from mne import io 
import os
from support_functions import get_file_prefixes, setup_logging, initialise_paths, clear_terminal
import gc
import logging
import json

paths = initialise_paths()
cleanpath = paths['cleanpath']
logspath = paths['logspath']

logger = setup_logging(logspath, 'clearn04')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
identify_twave_with_trapeziums = cfg["identify_twave_with_trapeziums"]
run_manual_ecg_insepction = cfg["run_manual_ecg_insepction"]

subjects = get_file_prefixes(cleanpath,'.fdt','.')

if not identify_twave_with_trapeziums:
    logger.info('Not performing identification of T-wave using trapeziums')
else:
    logger.info('Performing identification of T-wave using trapeziums')

if not run_manual_ecg_insepction:
    logger.info('Skipping manual ECG inspection')
else:
    logger.info('Performing manual ECG inspection')

gc.collect()
clear_terminal()
