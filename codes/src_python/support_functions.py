import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mne.viz import plot_topomap  
import pickle
import mne
from mne import combine_evoked
from scipy.signal import find_peaks
from itertools import product, combinations
import pandas as pd
import logging

def clear_terminal():
    os.system("cls" if os.name == "nt" else "clear")

def setup_logging(log_dir, script_name):
    """
    Configures logging for the current script.
    
    Args:
        log_dir (str): Directory where logs will be saved.
        script_name (str): Name of the script (used to name the log file).
    
    Returns:
        logging.Logger: Configured logger instance.
    """
    # Ensure the log directory exists
    os.makedirs(log_dir, exist_ok=True)
    
    # Define log file name
    log_file = os.path.join(log_dir, f"{script_name}_log.txt")
    
    # Configure logging
    logging.basicConfig(
        filename=log_file,
        filemode="w",
        level=logging.DEBUG,  # Log everything to the file
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Create a console handler for minimal output
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)  # Show only intermediate progress
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    logger = logging.getLogger()
    logger.addHandler(console_handler)

    return logger

def initialise_paths():
    """
    Initialise and ensure directories for commonly used paths.

    Returns
    -------
    paths : dict
        Dictionary containing initialized paths.
    """
    paths = {
        "root": os.path.join(os.pardir, os.pardir, "data"),
        "cleanpath": os.path.join(os.pardir, os.pardir, "data", "clean"),
        "preprocpath": os.path.join(os.pardir, os.pardir, "data", "preproc"),
        "eventspath": os.path.join(os.pardir, os.pardir, "data", "events"),
        "estimatespath": os.path.join(os.pardir, os.pardir, "data", "estimates"),
        "epochspath": os.path.join(os.pardir, os.pardir, "data", "epochs"),
        "GApath": os.path.join(os.pardir, os.pardir, "data", "GA"),
        "evokedspath": os.path.join(os.pardir, os.pardir, "data", "evokeds"),
        "tfrpath": os.path.join(os.pardir, os.pardir, "data", "tfr"),
        "pseudotrialpath": os.path.join(os.pardir, os.pardir, "data", "pseudotrials"),
        "derivativespath": os.path.join(os.pardir, os.pardir, "data", "derivatives"),
        "logspath": os.path.join(os.pardir, os.pardir, "data", "derivatives", "logs"),
        "vispath": os.path.join(os.pardir, os.pardir, "vis"),
        "cleaningvispath": os.path.join(os.pardir, os.pardir, "vis", "cleaning"),
        "latenciesvispath": os.path.join(os.pardir, os.pardir, "vis", "latencies"),
        "modelingvispath": os.path.join(os.pardir, os.pardir, "vis", "modeling"), # can drop?
        "pseudotrialvispath": os.path.join(os.pardir, os.pardir, "vis", "pseudotrials")
    }

    for path_name, path in paths.items():
        if not os.path.exists(path):
            os.makedirs(path)
            logging.info(f"Created directory: {path}")

    return paths

def get_file_prefixes(folderpath, endsWith, splitBy):
    
    # List all files in the data folder
    files = os.listdir(folderpath)
    
    # Extract the prefix (e.g., '01', '02', etc.) from each filename
    # Assumes filenames are in the format 'XX'+splitBy+endsWith where XX can vary
    prefixes = [filename.split(splitBy)[0] for filename in files if filename.endswith(endsWith)]
    
    # Sort the prefixes
    sorted_prefixes = sorted(prefixes)
    
    return sorted_prefixes


def find_next_same_cue_index(data, current_index):
    # Find the next trial with the same cue within the same block and subID
    
    current_subID = data.loc[current_index, 'subID']
    current_block = data.loc[current_index, 'block']
    current_cue = data.loc[current_index, 'trueAversive']
    
    # Start looking from the next row in the dataframe
    for next_index in range(current_index + 1, len(data)):
        if data.loc[next_index, 'subID'] == current_subID and data.loc[next_index, 'block'] == current_block:
            if data.loc[next_index, 'trueAversive'] == current_cue:
                return next_index
    return None

def extract_evoked_latency(evoked_data, evoked_type):
    # find frn, p3a, p3b peaks in evoked_data

    if evoked_type=='frn':
        
        cut_data = evoked_data.copy().pick(['FCz']).crop(0.2,0.37)
        latencies = cut_data.times[find_peaks(-cut_data.data[0],distance=100)[0]]
        
    elif evoked_type=='p3a':

        cut_data = evoked_data.copy().pick(['FCz']).crop(0.3,0.5)
        latencies = cut_data.times[find_peaks(cut_data.data[0],distance=100)[0]]

    elif evoked_type=='p3b':

        cut_data = evoked_data.copy().pick(['Pz']).crop(0.3,0.5)
        latencies = cut_data.times[find_peaks(cut_data.data[0],distance=100)[0]]

    del cut_data
    
    return latencies

def list_of_conditions(save_GA_main=False, save_individual_erps_main=False, 
                       save_GA_additional=False, save_individual_erps_additional=False,
                       save_GA_sys_dia=False, save_individual_erps_sys_dia=False, 
                       pseudotrial_conditions=False, compute_all_pseudotrials=False):

    rpeaks_main = ['changeIBI4_3', 'changeIBI5_4', 'changeIBI6_5', 'changeIBI7_6']
    rpeaks_control = ['changeIBI3_2', 'changeIBI2_1', 'changeIBI8_7']
    behaviors = ['outcome', 'absPE_2groups', 'expectation_3groups', 'learning_miniblock','choiceApproach']

    # for tfr and evokeds only:
    conditions = [['outcome'], ['sysDiaFeedb'], ['sysDiaFeedb', 'outcome']]
    conditions_ga = [['outcome'], ['sysDiaFeedb'], ['sysDiaFeedb', 'outcome']]
    conditions_ind = [['outcome'], ['sysDiaFeedb'], ['sysDiaFeedb', 'outcome']]
    conditions_pseudo = [['sysDiaFeedb']]
    
    # for tfr and evokeds only:
    if save_GA_main or save_individual_erps_main:
        conditions_hr_pe_main = [['sysDiaFeedb', rpeak+'_2groups'] for rpeak in rpeaks_main]
        conditions.extend(conditions_hr_pe_main)
        if save_individual_erps_main:
            conditions_ind.extend(conditions_hr_pe_main)        
        if save_GA_main:
            conditions_ga.extend(conditions_hr_pe_main)    
    # for pseudo
    if save_individual_erps_main or compute_all_pseudotrials:
        conditions_hr_pe_main = [['sysDiaFeedb', rpeak+'_2groups'] for rpeak in rpeaks_main]
        conditions_pseudo.extend(conditions_hr_pe_main)

    # for tfr and evokeds only:
    if save_GA_additional or save_individual_erps_additional:
        conditions_hr_pe_control1 = [['sysDiaFeedb', rpeak+'_2groups'] for rpeak in rpeaks_control]
        conditions.extend(conditions_hr_pe_control1)
        conditions_hr_pe_control2 = [['sysDiaFeedb', behavior, rpeak+'_2groups'] for rpeak in rpeaks_main for behavior in behaviors]
        conditions.extend(conditions_hr_pe_control2)
        if save_individual_erps_additional:
            conditions_ind.extend(conditions_hr_pe_control1)       
            conditions_ind.extend(conditions_hr_pe_control2)  
        if save_GA_additional:
            conditions_ga.extend(conditions_hr_pe_control1) 
            conditions_ga.extend(conditions_hr_pe_control2)
    # for pseudo
    if save_individual_erps_additional or compute_all_pseudotrials:
        conditions_hr_pe_control1 = [['sysDiaFeedb', rpeak+'_2groups'] for rpeak in rpeaks_control]
        conditions_pseudo.extend(conditions_hr_pe_control1)
        conditions_hr_pe_control2 = [['sysDiaFeedb', 'learning_miniblock', rpeak+'_2groups'] for rpeak in rpeaks_main]
        conditions_pseudo.extend(conditions_hr_pe_control2)
        conditions_hr_pe_control3 = [[rpeak+'_2groups'] for rpeak in rpeaks_main]
        conditions_pseudo.extend(conditions_hr_pe_control3)
        #conditions_pseudo.extend([['sysDiaFeedb','beforeFeedbAveIBI_2groups'],['sysDiaFeedb','afterFeedbAveIBI_2groups']])

    # for tfr and evokeds only:
    if save_GA_sys_dia or save_individual_erps_sys_dia:
        conditions_sys_dia = [['sysDiaFeedb', 'outcome', 'learning_miniblock']]
        conditions.extend(conditions_sys_dia)
        if save_individual_erps_sys_dia:
            conditions_ind.extend(conditions_sys_dia)  
        if save_GA_sys_dia:
            conditions_ga.extend(conditions_sys_dia)
    # for pseudo
    if save_individual_erps_sys_dia or compute_all_pseudotrials:
        conditions_sys_dia = [['sysDiaFeedb','learning_miniblock'],['learning_miniblock']]
        conditions_pseudo.extend(conditions_sys_dia)

    if pseudotrial_conditions:
        conditions = conditions_pseudo
    else:
        conditions = [conditions, conditions_ga, conditions_ind]
    
    return conditions

def get_miniblock_intervals(feedb_latencies, time_range_blocks, sfreq):

    learning_miniblock = np.tile(np.repeat(['first', 'next', 'reversal', 'last'],25),3)
    
    miniblock_intervals = []
    for block_idx, time_range in enumerate(time_range_blocks):
        # Filter feedback latencies and miniblocks for the current block
        block_latencies = feedb_latencies[
            (feedb_latencies > time_range[0]) & (feedb_latencies < time_range[1])
        ]
        block_labels = learning_miniblock[:len(block_latencies)] 
    
        # Define miniblock start and end points
        block_boundaries = []
        current_start = time_range[0]  # Start of the block
        for label in ['first', 'next', 'reversal', 'last']:
            last_latency = block_latencies[np.where(block_labels == label)[0][-1]]+sfreq*2 # trial: feedback+2seconds
            block_boundaries.append((current_start, last_latency, label))
            current_start = last_latency+1 # Update for the next miniblock
    
        # 'last' miniblock ends at the block end
        block_boundaries[-1] = (block_boundaries[-1][0], time_range[1], 'last')
        miniblock_intervals.extend(block_boundaries)

    return miniblock_intervals

def get_ibi_group_thresholds(split_IBI_into_groups_for_pseudotrials, thresholdpath = 'thresholdpath'):

    if split_IBI_into_groups_for_pseudotrials == 'all_ind_ibi':  
        _, thresh_ibi  = pd.qcut(IBIs, 2, retbins=True)
        _, thresh_change_ibi = pd.qcut(change_IBI, 2, retbins=True)
    
    elif split_IBI_into_groups_for_pseudotrials == 'slow_speed':   
        thresh_ibi  = 2
        thresh_change_ibi = np.array([-10000, 0, 10000])
    
    elif split_IBI_into_groups_for_pseudotrials == 'thresholds_from_data':   
        # load pickle that contains a dictionary of thresholds per subject, keys - variables
        # or just load df and rerun for all needed vars
        thresh_ibi  = thresholdpath
        thresh_change_ibi = thresholdpath
        
    else: # by default 'all_ind_ibi_within_iteration'    
        thresh_ibi  = 2
        thresh_change_ibi = 2

    return thresh_ibi, thresh_change_ibi

def create_metadata_for_pseudotrials(new_feedback_latencies, Rpeak_latencies, IBIs, change_IBI, miniblock_intervals, thresh_ibi, thresh_change_ibi):

    ibi_labels = ['IBI' + str(i) for i in range(1, 9)]
    change_labels = ['change' + ibi_labels[i] + '_' + str(i) for i in range(1,len(ibi_labels))]
    change_group_labels = ['change' + ibi_labels[i] + '_' + str(i) + '_2groups' for i in range(1,len(ibi_labels))]
    cols = ibi_labels + change_labels + change_group_labels + ['sysDiaFeedb', 'phase', 'learning_miniblock', 'beforeFeedbAveIBI', 'afterFeedbAveIBI', 'beforeFeedbAveIBI_2groups', 'afterFeedbAveIBI_2groups']
    meta = pd.DataFrame(columns=cols, index=range(300))

    for feedb in range(len(new_feedback_latencies)):    
        
        ibi_before_feedb = IBIs[Rpeak_latencies<new_feedback_latencies[feedb]][-4:]
        ibi_after_feedb = IBIs[Rpeak_latencies>new_feedback_latencies[feedb]][:4]
        change_ibi_before_feedb = change_IBI[Rpeak_latencies<new_feedback_latencies[feedb]][-3:]
        change_ibi_after_feedb = change_IBI[Rpeak_latencies>new_feedback_latencies[feedb]][:4]
        
        meta.loc[feedb, ibi_labels] = np.concatenate([ibi_before_feedb, ibi_after_feedb])
        meta.loc[feedb, change_labels] = np.concatenate([change_ibi_before_feedb, change_ibi_after_feedb])
        meta.loc[feedb, 'beforeFeedbAveIBI'] = np.mean(ibi_before_feedb)
        meta.loc[feedb, 'afterFeedbAveIBI'] = np.mean(ibi_after_feedb)
            
        sampleDistToFeedb = new_feedback_latencies[feedb] - Rpeak_latencies[Rpeak_latencies<new_feedback_latencies[feedb]][-1]
        meta.loc[feedb, 'phase'] = (sampleDistToFeedb/ibi_after_feedb[0])*360
    
    meta.loc[meta.phase<120,'sysDiaFeedb'] = 'Systole'
    meta.loc[meta.phase>240,'sysDiaFeedb'] = 'Diastole'

    pseudotrial_labels = []
    for pseudotrial in new_feedback_latencies:
        for start, end, label in miniblock_intervals:
            if start <= pseudotrial <= end:
                pseudotrial_labels.append(label)
                break

    meta['learning_miniblock'] = pseudotrial_labels

    create_heart_groups = change_labels + ['beforeFeedbAveIBI', 'afterFeedbAveIBI']
    
    for var in create_heart_groups:
    
        if var in change_labels:
            thresh = thresh_change_ibi
        else:
            thresh = thresh_ibi
    
        if isinstance(thresh, np.ndarray):
            meta[var+'_2groups'] = pd.cut(meta[var], bins=thresh, labels=['low','high']).values
            print('for variable %s, the thresholds were already provided. those are %s'%(var,thresh))
        elif isinstance(thresh, dict):
            meta[var+'_2groups'] = pd.cut(meta[var], bins=thresh[var], labels=['low','high']).values
            print('for variable %s, the thresholds were provided in a dict. those are %s'%(var,thresh[var]))
        else: # thresh ==2:
            meta[var+'_2groups'] = pd.qcut(meta[var], 2, labels=['low','high']).values
            print('for variable %s, the thresholds were computed based on current iteration. there are 2 bins'%(var))

    return meta

def get_event_dict(new_feedback_latencies):
    
    new_event_dict = {'feedback':1}
    new_triggers = np.hstack((np.repeat(new_event_dict['feedback'],len(new_feedback_latencies)))) 
    zeros = np.repeat(0,len(new_triggers))
    new_events = np.vstack((new_feedback_latencies,zeros,new_triggers)).T
    new_events = new_events[new_events[:, 0].argsort()]

    return new_events, new_event_dict

def define_evokeds(epochs, sub_df, variables, bad_epochs=False, baseline=False, baseline_window=None, baseline_mode=None):
    # epochs - epochs of this subject
    # sub_df - metadata, same length as epochs
    # variables - list of variables (array of strings) to use for defining evokeds
    # bad_epochs - bool, drop bad epochs before averaging
    # baseline - bool, apply baseline before averaging
    # baseline window - tuple, in seconds

    if baseline:
        if isinstance(epochs, mne.time_frequency.EpochsTFR):
            epochs.apply_baseline(baseline_window, mode=baseline_mode)
        else:
            epochs.apply_baseline(baseline_window)
                
    # for each variable, determine factors
    dict_levels = {} 
    for variable in variables:
        dict_levels[variable] = sub_df[variable].astype('category').cat.categories.tolist()

    variable_levels = list(dict_levels.values())

    all_combinations = product(*variable_levels)
    dict_evokeds = {}
    dict_num_ave = {}
    
    for combination in all_combinations:
        
        interaction_string = "_".join([f"{var}_{level}" for var, level in zip(variables, combination)])
           
        mask = pd.Series(True, index=sub_df.index)
        for var, level in zip(variables, combination):
            mask &= (sub_df[var] == level)

        if bad_epochs:
             mask &= (sub_df['bad_epoch']!=1)
        
        dict_num_ave[interaction_string] = sum(mask)
        
        try:
            if isinstance(epochs, mne.time_frequency.EpochsTFR):
                dict_evokeds[interaction_string] = epochs[mask].average()
            else:
                dict_evokeds[interaction_string] = epochs[mask].average(picks='all')
        except:
            print(f"Error for interaction: {interaction_string}")
            dict_evokeds[interaction_string] = None

    return dict_evokeds, dict_num_ave

def calculate_evoked_differences(dict_evokeds, variables, subtraction_order=None):
    """
    Calculate evoked response differences based on existing keys in dict_evokeds,
    using a specified subtraction order if provided. Prioritizes subtracting from 'high' or 'first'
    if those levels are present; otherwise, it uses sorted order.

    Parameters:
        - dict_evokeds: dictionary of evoked responses for all conditions
        - variables: list of variables to consider for interactions
        - subtraction_order: dictionary specifying level pairs for each variable for subtraction,
                             e.g., {"levels": [["low", "high"], ["low", "medium"]], "sysDiaFeedb": [["Diastole", "Systole"]]}
    
    Returns:
        - dict_differences: dictionary with evoked differences
    """
    dict_differences = {}
    
    # Loop over each variable and calculate differences
    for variable in variables:
        # Identify levels of the current variable by scanning the keys in dict_evokeds
        levels = sorted(set(key.split(f"{variable}_")[1].split('_')[0] for key in dict_evokeds.keys() if f"{variable}_" in key))
        
        # Generate all possible pairs of levels
        pairs = list(combinations(levels, 2))

        # Adjust pairs based on `subtraction_order`, `high`, or `first`
        adjusted_pairs = []
        for first_level, second_level in pairs:
            # Check if specific order is defined in subtraction_order
            if subtraction_order and variable in subtraction_order:
                if [first_level, second_level] in subtraction_order[variable]:
                    # Keep order as is
                    adjusted_pairs.append((first_level, second_level))
                elif [second_level, first_level] in subtraction_order[variable]:
                    # Reverse order
                    adjusted_pairs.append((second_level, first_level))
                else:
                    # If not specified, use default ordering
                    if "high" in (first_level, second_level):
                        adjusted_pairs.append((second_level, first_level) if first_level == "high" else (first_level, second_level))
                    elif "first" in (first_level, second_level):
                        adjusted_pairs.append((second_level, first_level) if first_level == "first" else (first_level, second_level))
                    else:
                        # Default sorted order
                        adjusted_pairs.append((first_level, second_level))
            else:
                # Apply default rules if subtraction_order not provided
                if "high" in (first_level, second_level):
                    adjusted_pairs.append((second_level, first_level) if first_level == "high" else (first_level, second_level))
                elif "first" in (first_level, second_level):
                    adjusted_pairs.append((second_level, first_level) if first_level == "first" else (first_level, second_level))
                else:
                    # Default sorted order
                    adjusted_pairs.append((first_level, second_level))

        # Calculate differences for each adjusted pair
        for first_level, last_level in adjusted_pairs:
            # Now, find pairs of keys that differ only by the current variable's level
            for key in dict_evokeds.keys():
                if f"{variable}_{first_level}" in key:
                    # Construct the corresponding key for the last level by replacing the level in the key
                    counterpart_key = key.replace(f"{variable}_{first_level}", f"{variable}_{last_level}")
                    
                    # Check that both keys exist in dict_evokeds and are not None
                    evoked_first = dict_evokeds.get(key)
                    evoked_last = dict_evokeds.get(counterpart_key)

                    if evoked_first is not None and evoked_last is not None:
                        
                        if key.startswith("tfr_") and counterpart_key.startswith("tfr_"):
                            base_key = key[4:]  # Remove 'tfr_'
                            base_counterpart_key = counterpart_key[4:]  # Remove 'tfr_'
                            diff_key = f"tfr_dif_{base_counterpart_key}_{base_key}"
                        else:
                            diff_key = f"dif_{counterpart_key}_{key}"

                        # Check if the difference has already been computed
                        if diff_key in dict_differences:
                            logging.info(f"Skipping duplicate computation for {diff_key}")
                            continue

                        # Compute the difference (last - first)
                        dict_differences[diff_key] = combine_evoked([evoked_last, evoked_first], [1, -1])
                        logging.info(f"Computed difference for: {diff_key}")
                    else:
                        # Log if one of the evoked responses is missing
                        if evoked_first is None:
                            logging.warning(f"Missing data for {key}")
                        if evoked_last is None:
                            logging.warning(f"Missing data for {counterpart_key}")

    return dict_differences

def load_evoked_response(subject, key, files, pseudotrialpath):
    """
    Load a single evoked response for a given subject and key.

    Args:
        subject (str): Subject ID.
        key (str): The key (condition) to process.
        files (list): List of available files.
        pseudotrialpath (str): Path to the pseudotrial directory.

    Returns:
        evoked: The loaded evoked response, or None if loading fails.
    """
    for ext in ["_ave.hdf5", "_ave.fif"]:
        filename = f"{subject}_{key}{ext}"
        if filename in files:
            file_path = os.path.join(pseudotrialpath, filename)
            try:
                if ext == "_ave.hdf5":
                    evoked = mne.time_frequency.read_tfrs(file_path)
                else:
                    evoked = mne.read_evokeds(file_path)[0]
                logging.info(f"Loaded evoked response for {filename}")
                return evoked
            except Exception as e:
                logging.error(f"Failed to load {filename}: {e}")
                return None
    logging.warning(f"No matching file found for {subject}_{key}")
    return None


def plot_significant_clusters(good_cluster_inds, clusters, T_obs, evokeds_dif, times, ch_info, vispath, conditions):

    """
    Plots significant clusters and saves each cluster plot as a PNG.

    Parameters
    ----------
    good_cluster_inds : list
        Indices of significant clusters.
    clusters : list
        List of clusters from the statistical test.
    T_obs : ndarray
        Observed T-statistics for each time point and channel.
    grand_ave : ndarray
        Grand average data across subjects.
    times : ndarray
        Time points corresponding to the data.
    evokeds_win : list of Evoked
        Evoked data used for plotting topomaps.
    vispath : str
        Directory path where to save the plots.
    """
    
    grand_ave = evokeds_dif.mean(axis=0)
    
    for i_clu, clu_idx in enumerate(good_cluster_inds):
        # Unpack cluster information, get unique indices
        time_inds, space_inds = np.squeeze(clusters[clu_idx])
        ch_inds = np.unique(space_inds)
        time_inds = np.unique(time_inds)

        # Get topography for F stat
        f_map = T_obs[time_inds, ...].mean(axis=0)

        # Get signals at significant sensors
        signals = grand_ave[..., ch_inds].mean(axis=1)
        sig_times = times[time_inds]

        # Create spatial mask
        mask = np.zeros((f_map.shape[0], 1), dtype=bool)
        mask[ch_inds, :] = True

        # Initialize figure
        fig, ax_topo = plt.subplots(1, 1, figsize=(10, 3))
        title = 'Cluster #{0}'.format(i_clu + 1)
        fig.suptitle(title, fontsize=14)

        # Plot average test statistic and mark significant sensors
        image, _ = plot_topomap(f_map, ch_info, mask=mask, axes=ax_topo,
                                cmap='Reds', show=False, sphere=(0.00, 0.02, 0.00, 0.09))

        # Advanced matplotlib for showing image with figure and colorbar
        divider = make_axes_locatable(ax_topo)

        # Add axes for colorbar
        ax_colorbar = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(image, cax=ax_colorbar)
        ax_topo.set_xlabel('Averaged T-map ({:0.1f} - {:0.1f} ms)'.format(
            *sig_times[[0, -1]]
        ))

        # Add new axis for time courses and plot time courses
        ax_signals = divider.append_axes('right', size='300%', pad=1.2)
        ax_signals.plot(times, signals)

        # Add information
        ax_signals.axvline(0, color='k', linestyle=':')
        ax_signals.set_xlim([times[0], times[-1]])
        ax_signals.set_xlabel('time [ms]')
        ax_signals.set_ylabel('effect size [V]')

        # Plot significant time range
        ymin, ymax = ax_signals.get_ylim()
        ax_signals.fill_betweenx((ymin, ymax), sig_times[0], sig_times[-1],
                                 color='orange', alpha=0.3)
        ax_signals.set_ylim(ymin, ymax)

        # Save the figure instead of displaying it
        fig_name = f'cluster_{i_clu + 1}'+ '_' + conditions[0] + '_' + conditions[1] + '.png'
        fig_path = os.path.join(vispath, fig_name)
        fig.savefig(fig_path)
        plt.close(fig)  # Close the figure after saving to avoid display

def prepare_for_clustering(conditions, evokedspath):   
    filepath1 = os.path.join(evokedspath, conditions[0]+'.pkl')

    # Load the pickled evokeds list
    with open(filepath1, 'rb') as f:
        evokeds_1 = pickle.load(f)
    
    ch_info = evokeds_1[0].copy().pick(['eeg']).info
    times = evokeds_1[0].times
    
    filepath2 = os.path.join(evokedspath, conditions[1]+'.pkl')
        
    # Load the pickled evokeds list
    with open(filepath2, 'rb') as f:
        evokeds_2 = pickle.load(f)
    
    evokeds_dif = []
    
    for sub in range(len(evokeds_1)):
    
        evokeds_dif.append(mne.combine_evoked([evokeds_1[sub],evokeds_2[sub]],[1,-1]))
        evokeds_dif[-1].pick(['eeg'])
        evokeds_dif[-1] = evokeds_dif[-1].data

    return evokeds_dif, ch_info, times

def preprocess_pseudotrial_filename(filename, is_tfr=False):
    """Remove subject prefix and file suffix from pseudotrial filename."""
    if is_tfr:
        return "_".join(filename.split("_")[2:]).replace("_ave.hdf5", "")
    return "_".join(filename.split("_")[1:]).replace("_ave.fif", "")

def find_best_match(evoked_key, pseudotrial_keys):
    """
    Find the pseudotrial key(s) that match the most keywords from evoked_key.

    Parameters:
    ----------
    evoked_key : str
        The evoked condition key.
    pseudotrial_keys : list of str
        List of available pseudotrial keys.

    Returns:
    -------
    list
        A list of all best matches.
    """
    evoked_keywords = evoked_key.split("_")
    best_matches = []  # Track multiple best matches
    best_score = 0

    for pseudotrial_key in pseudotrial_keys:
        pseudotrial_keywords = pseudotrial_key.split("_")

        # Exclude pseudotrials containing keywords not in evoked_key
        if any(keyword not in evoked_keywords for keyword in pseudotrial_keywords):
            continue

        # Count the number of matching keywords
        match_score = sum(1 for kw in pseudotrial_keywords if kw in evoked_keywords)

        if match_score > best_score:
            best_score = match_score
            best_matches = [pseudotrial_key]
        elif match_score == best_score:
            best_matches.append(pseudotrial_key)

    return best_matches

def extract_neural_data(df, pick_sub, epochs, channels, roi, time_roi, freq=None, freq_band=None, dont_rescale=None):

    if isinstance(epochs, mne.time_frequency.EpochsTFR): 
        for channel in channels:
            if dont_rescale==True:
                df.loc[pick_sub, channel+'__'+freq+'__'+roi] = epochs.get_data(picks=channel, fmin=freq_band[0], fmax=freq_band[1], tmin=time_roi[0], tmax=time_roi[1])[:,0,:,:].mean(axis=1).mean(axis=1)
            else:
                df.loc[pick_sub, channel+'__'+freq+'__'+roi] = epochs.get_data(picks=channel, fmin=freq_band[0], fmax=freq_band[1], tmin=time_roi[0], tmax=time_roi[1])[:,0,:,:].mean(axis=1).mean(axis=1)*1e12 # converted to mV2 (it was in V2)
    elif isinstance(epochs, mne.Evoked):
        for channel in channels:
            if dont_rescale==True:
                df.loc[pick_sub, channel+'__'+roi] = epochs.get_data(picks=channel, tmin=time_roi[0], tmax=time_roi[1])[0,:].mean()
            else:
                df.loc[pick_sub, channel+'__'+roi] = epochs.get_data(picks=channel, tmin=time_roi[0], tmax=time_roi[1])[0,:].mean()*1e6 # converted to mV (it was in V)            
    elif isinstance(epochs, mne.time_frequency.AverageTFR):
        for channel in channels:
            if dont_rescale==True:
                df.loc[pick_sub, channel+'__'+freq+'__'+roi] = epochs.get_data(picks=channel, fmin=freq_band[0], fmax=freq_band[1], tmin=time_roi[0], tmax=time_roi[1])[0,:,:].mean()
            else:
                df.loc[pick_sub, channel+'__'+freq+'__'+roi] = epochs.get_data(picks=channel, fmin=freq_band[0], fmax=freq_band[1], tmin=time_roi[0], tmax=time_roi[1])[0,:,:].mean()*1e12 # converted to mV2 (it was in V2)
    else:
        for channel in channels:
            if dont_rescale==True:
                df.loc[pick_sub, channel+'__'+roi] = epochs.get_data(picks=channel, tmin=time_roi[0], tmax=time_roi[1], copy=True)[:,0,:].mean(axis=1)
            else:
                df.loc[pick_sub, channel+'__'+roi] = epochs.get_data(picks=channel, tmin=time_roi[0], tmax=time_roi[1], copy=True)[:,0,:].mean(axis=1)*1e6 # converted to mV (it was in V)

def extract_neural_data_multiple_roi(df, pick_sub, epochs, ch_ecg, ch_eog, ch_resp, ch_of_interest, eeg_of_interest, save_ecg, save_eog, save_breathing, rois, rename_roi="", baseline_window=None, baseline_window_tfr=None, freqs_of_interest=None, dont_rescale=None):

    renamed_rois = {f"{key}{rename_roi}": value for key, value in rois.items()}
    
    if isinstance(epochs, mne.time_frequency.EpochsTFR) or isinstance(epochs, mne.time_frequency.AverageTFR): 

        for freq, freq_band in freqs_of_interest.items():
        
            for roi, time_roi in renamed_rois.items():           
                if roi.startswith('p3b'):
                    channel = ['Pz']
                else:
                    channel = ['FCz']
                    
                extract_neural_data(df, pick_sub, epochs, channel, roi, time_roi, freq, freq_band, dont_rescale)

            if baseline_window:
                extract_neural_data(df, pick_sub, epochs, eeg_of_interest, 'baseline'+rename_roi, baseline_window, freq, freq_band, dont_rescale)
            if baseline_window_tfr:
                extract_neural_data(df, pick_sub, epochs, eeg_of_interest, 'baseline_tfr'+rename_roi, baseline_window_tfr, freq, freq_band, dont_rescale)
            
    else:
    
        for roi, time_roi in renamed_rois.items():            
            if roi.startswith('p3b'):
                channel = ['Pz']
            else:
                channel = ['FCz']
    
            extract_neural_data(df, pick_sub, epochs, channel, roi, time_roi, freq=None, freq_band=None, dont_rescale=dont_rescale)
            if save_ecg:
                extract_neural_data(df, pick_sub, epochs, ch_ecg, roi, time_roi, freq=None, freq_band=None, dont_rescale=dont_rescale)
            if save_eog:
                extract_neural_data(df, pick_sub, epochs, ch_eog, roi, time_roi, freq=None, freq_band=None, dont_rescale=dont_rescale)
            if save_breathing:
                extract_neural_data(df, pick_sub, epochs, ch_resp, roi, time_roi, freq=None, freq_band=None, dont_rescale=dont_rescale)
    
        if baseline_window:
            extract_neural_data(df, pick_sub, epochs, ch_of_interest, 'baseline'+rename_roi, baseline_window, freq=None, freq_band=None, dont_rescale=dont_rescale)
        if baseline_window_tfr:
            extract_neural_data(df, pick_sub, epochs, ch_of_interest, 'baseline_tfr'+rename_roi, baseline_window_tfr, freq=None, freq_band=None, dont_rescale=dont_rescale)

def process_files(file_path, df_num_ave, pseudotrialpath, subjects, gapath, is_tfr=False):
    
    pseudotrial_files = []
    if is_tfr:
        pseudotrial_files = list(set(preprocess_pseudotrial_filename(f, is_tfr=True) for f in os.listdir(pseudotrialpath) if f.endswith("_ave.hdf5")))
    else:
        pseudotrial_files = list(set(preprocess_pseudotrial_filename(f) for f in os.listdir(pseudotrialpath) if f.endswith("_ave.fif")))

    file_list = [f for f in os.listdir(file_path) if f.endswith(".pkl") and not f.startswith(("dif_", "baseline_", "long_")) and "_corrected" not in f]

    for file in file_list:
        condition_key = file.replace(".pkl", "")

        if condition_key not in df_num_ave.columns:
            logging.warning(f"Condition {condition_key} not found in df_num_ave columns")
            continue

        best_match_keys = find_best_match(condition_key, pseudotrial_files)

        if not best_match_keys:
            logging.warning(f"No suitable pseudotrial found for {condition_key}")
            continue

        if len(best_match_keys)>1:
            logging.warning(f"Several equally suitable pseudotrials found for {condition_key}: {best_match_keys}")

        valid_subjects = [subject_prefix for subject_prefix in subjects if df_num_ave.at[int(subject_prefix), condition_key] > 0]

        # Process valid subjects
        full_path = os.path.join(file_path, file)
        with open(full_path, 'rb') as f:
            data_list = pickle.load(f)  

        for best_match_key in best_match_keys:
            corrected_list = []
            for subject_idx, subject_prefix in enumerate(valid_subjects):
                if subject_idx >= len(data_list):
                    logging.warning(f"Subject index {subject_idx} out of range for data_list. Skipping subject {subject_prefix}.")
                    continue
    
                data = data_list[subject_idx]
    
                pseudotrial_filename = f"{subject_prefix}{'_tfr_' if is_tfr else '_'}{best_match_key}_ave.{'hdf5' if is_tfr else 'fif'}"
                pseudotrial_path = os.path.join(pseudotrialpath, pseudotrial_filename)
    
                if not os.path.exists(pseudotrial_path):
                    logging.warning(f"Pseudotrial file {pseudotrial_filename} not found for subject {subject_prefix}")
                    continue
    
                pseudotrial = (mne.time_frequency.read_tfrs(pseudotrial_path) if is_tfr else mne.read_evokeds(pseudotrial_path)[0])
    
                if len(data.times) != len(pseudotrial.times):
                    data.crop(tmin=pseudotrial.times[0]+0.00000001, tmax=pseudotrial.times[-1]-0.0000001)
                    logging.info('Times in initial data were imprecise. Cropping.')
    
                # Perform the subtraction
                corrected_data = mne.combine_evoked([data, pseudotrial], weights=[1, -1])
                corrected_list.append(corrected_data)
    
                # baseline correct pseudo-corrected evoked here if needed (in this case, also need separate corrected_list_baseline to fill out)
                # not 100% correct: better to baseline-correct pseudotrials already and subtract from baseline-corrected individual evokeds?
    
            # Save the corrected list as a pickle
            corrected_pickle_path = os.path.join(file_path, f"{condition_key}_corrected_{best_match_key}.pkl")
            with open(corrected_pickle_path, 'wb') as f:
                pickle.dump(corrected_list, f)
    
            logging.info(f"Saved corrected {'TFR' if is_tfr else 'Evoked'} list for condition {condition_key} with pseudotrial {best_match_key} at {corrected_pickle_path}")
    
            try:
                save_group_erp(f"{condition_key}_corrected_{best_match_key}", corrected_list, gapath, save_ga=True, save_csv_with_ci=True, channels=['FCz', 'Pz'], return_ga=False, is_tfr=is_tfr)
            except Exception as e:
                logging.error(f"Error saving grand average for corrected {condition_key} with pseudotrial {best_match_key}: {e}")
    
            clear_output(wait=True)
            # clear_terminal()

def compute_condition_differences(file_path, df_num_ave, subjects, gapath, is_tfr=False):
    # Get lists of corrected and difference files
    corrected_files = [f for f in os.listdir(file_path) if f.endswith(".pkl") and "_corrected_" in f]
    difference_files = [f for f in os.listdir(file_path) if f.startswith("dif_") and f.endswith(".pkl") and "_corrected" not in f]

    processed_pairs = set()

    for dif_file in difference_files:
        # Extract condition keys from the difference file name
        file_base = dif_file.replace("dif_", "").replace(".pkl", "")
        matched = False
        for corrected_file in corrected_files:
            # Extract the condition and pseudotrial from the corrected file
            condition_key1, pseudotrial_key1 = corrected_file.replace(".pkl", "").split("_corrected_")
            # check that condition_key1 is in file_base
            condition_key2_candidate = file_base[len(condition_key1) + 1:]
            first_part1 = condition_key1.split("_")[0]
            first_part2 = condition_key2_candidate.split("_")[0]
            if first_part1 == first_part2:
                condition_key2 = condition_key2_candidate
                matched = True
                break
        if not matched:
            logging.warning(f"Could not resolve condition pair for {dif_file}")
            continue

        if (condition_key1, condition_key2) in processed_pairs or (condition_key2, condition_key1) in processed_pairs:
            continue

        # Filter corrected files for the second condition
        condition2_corrected_files = [f for f in corrected_files if f.startswith(f"{condition_key2}_corrected_")]

        for corrected_file2 in condition2_corrected_files:
            _, pseudotrial_key2 = corrected_file2.replace(".pkl", "").split("_corrected_")

            corrected_file1 = f"{condition_key1}_corrected_{pseudotrial_key1}.pkl"
            corrected_file2 = f"{condition_key2}_corrected_{pseudotrial_key2}.pkl"

            if corrected_file1 not in corrected_files or corrected_file2 not in corrected_files:
                logging.warning(f"Missing corrected files for conditions {condition_key1} or {condition_key2}")
                continue

            # Check pseudotrial differences
            if pseudotrial_key1 == pseudotrial_key2:
                logging.info(f"Skipping computation: pseudotrial keys are identical ({pseudotrial_key1})")
                continue

            if len(set(pseudotrial_key1.split("_")).symmetric_difference(set(pseudotrial_key2.split("_")))) > 2:
                logging.info(f"Not computing differences for {condition_key1} and {condition_key2} with pseudo keys {pseudotrial_key1} and {pseudotrial_key2} - they are too different")
                continue

            # Load the corrected files
            with open(os.path.join(file_path, corrected_file1), 'rb') as f:
                list1 = pickle.load(f)
            with open(os.path.join(file_path, corrected_file2), 'rb') as f:
                list2 = pickle.load(f)

            # Identify valid subjects for both conditions
            valid_subjects1 = [subject_prefix for subject_prefix in subjects if df_num_ave.at[int(subject_prefix), condition_key1] > 0]
            valid_subjects2 = [subject_prefix for subject_prefix in subjects if df_num_ave.at[int(subject_prefix), condition_key2] > 0]
            common_subjects = list(set(valid_subjects1).intersection(valid_subjects2))

            # Filter lists to only include common subjects
            filtered_list1 = [list1[valid_subjects1.index(subject)] for subject in common_subjects]
            filtered_list2 = [list2[valid_subjects2.index(subject)] for subject in common_subjects]

            if len(filtered_list1) != len(filtered_list2):
                logging.warning(f"Filtered mismatch in number of subjects for {condition_key1} and {condition_key2}")
                continue

            # Compute differences for each subject
            difference_list = []
            for idx, (data1, data2) in enumerate(zip(filtered_list1, filtered_list2)):
                diff = mne.combine_evoked([data1, data2], weights=[1, -1])
                difference_list.append(diff)

            # Save the difference list
            dif_corrected_file = f"dif_{condition_key1}_corrected_{pseudotrial_key1}__{condition_key2}_corrected_{pseudotrial_key2}.pkl"
            with open(os.path.join(file_path, dif_corrected_file), 'wb') as f:
                pickle.dump(difference_list, f)

            logging.info(f"Saved corrected differences for {condition_key1} and {condition_key2} at {dif_corrected_file}")

            try:
                save_group_erp(f"dif_{condition_key1}_corrected_{pseudotrial_key1}__{condition_key2}_corrected_{pseudotrial_key2}",
                               difference_list, gapath, save_ga=True, save_csv_with_ci=True, channels=['FCz', 'Pz'], return_ga=False, is_tfr=is_tfr)
            except Exception as e:
                logging.error(f"Error saving grand average for corrected {condition_key1} and {condition_key2}: {e}")

            clear_output(wait=True)
          
        processed_pairs.add((condition_key1, condition_key2))

def save_channels_to_csv(key, GApath, grand_avg, channels, evoked_list=None):
    
    data_dict = {"Time (s)": np.round(grand_avg.times, decimals=3)}

    for channel in channels:
        if channel not in grand_avg.ch_names:
            logging.warning(f"Channel '{channel}' not found in grand average. Skipping.")
            continue

        if evoked_list:
            
            # Stack data from individual evoked responses for the current channel
            channel_data = np.array([evoked.data[evoked.ch_names.index(channel)] for evoked in evoked_list])

            # Compute mean and 95% CI (based on standard error)
            def compute_ci(data):
                mean = data.mean(axis=0)
                sem = data.std(axis=0, ddof=1) / np.sqrt(data.shape[0])
                ci_upper = mean + 1.96 * sem
                ci_lower = mean - 1.96 * sem
                return mean, ci_lower, ci_upper

            mean, ci_lower, ci_upper = compute_ci(channel_data)

            # Convert to microvolts (from volts)
            ci_lower *= 1e6
            ci_upper *= 1e6

            # Add data to dictionary
            data_dict[f"{channel}_CI_lower (uV)"] = ci_lower
            data_dict[f"{channel}_CI_upper (uV)"] = ci_upper

        else:
            
            mean = grand_avg.get_data(picks=channel)[0]
            
        mean *= 1e6
        data_dict[f"{channel}_mean (uV)"] = mean

    # Prepare DataFrame for CSV
    df = pd.DataFrame(data_dict)

    # Save CSV     
    csv_path = os.path.join(GApath, f"{key}_Channels{'_CI' if evoked_list else ''}.csv")
    df.to_csv(csv_path, index=False)
    logging.info(f"Saved data {'with 95% CI ' if evoked_list else ''}for channels {channels} for '{key}' at: {csv_path}")


def save_group_erp(key, evoked_list, GApath, save_ga=True, save_csv_with_ci=True, channels=['FCz', 'Pz'], return_ga=False, is_tfr=False):
    """
    Save grand average (GA) and optionally extract signals from specified channels with their 95% CI.

    Parameters
    ----------
    key : str
        Condition key to label the output files.
    evoked_list : list of mne.Evoked
        List of Evoked objects.
    GApath : str
        Path to save the outputs.
    save_ga : bool
        If True, saves the grand average as a .fif or .hdf5 file.
    save_csv_with_ci : bool
        If True, saves specified channel signals along with 95% confidence intervals as a CSV file.
    channels : list of str
        List of channel names to extract data from.
    return_ga : bool
        If True, returns the computed grand average.
    is_tfr : bool
        If True, handles time-frequency representations and adjusts output file formats.
    """
    
    evoked_list = [x for x in evoked_list if x is not None]
    
    grand_avg = mne.grand_average(evoked_list)

    if save_ga:
        ga_path = os.path.join(GApath, f"{'tfr_' if is_tfr else ''}{key}_GA.{'hdf5' if is_tfr else 'fif'}")
        grand_avg.save(ga_path, overwrite=True, verbose=False)
        logging.info(f"Saved GA for condition '{key}' at: {ga_path}")

    if is_tfr:
        save_csv_with_ci = False

    if save_csv_with_ci:

        save_channels_to_csv(key, GApath, grand_avg, channels, evoked_list)

    if return_ga:
        return grand_avg
    else:
        del grand_avg
