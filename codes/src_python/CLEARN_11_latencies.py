import numpy as np
import pandas as pd
import os
import mne
from support_functions import extract_evoked_latency, setup_logging, initialise_paths, clear_terminal
import pickle
from collections import defaultdict
from scipy.stats import ttest_rel
from matplotlib import pyplot as plt
import gc
import logging
import json

paths = initialise_paths()
root = paths['root']
evokedspath = paths['evokedspath']
GApath = paths['GApath']
derivativespath = paths['derivativespath']
logspath = paths['logspath']
vispath = paths['vispath']

logger = setup_logging(logspath, 'clearn11')

config_path = os.path.join(os.pardir,'global_variables.json')
with open(config_path,'r') as f:
    cfg = json.load(f)
significance_level_for_peak_differences = 0.1
window_length = 0.05
check_valence_phase_interactions = True
save_vis = True
round_num = 3

erp_types = ['frn','p3a','p3b'] 
peak_latencies = defaultdict(list)
final_latencies = defaultdict(list)

lat_win = defaultdict(list)
lat_lose = defaultdict(list)
GA = defaultdict(list)

GA['overall'] = mne.read_evokeds(os.path.join(GApath, 'clean_GA.fif'))[0]
GA['Win'] = mne.read_evokeds(os.path.join(GApath, 'outcome_Win_GA.fif'))[0]
GA['Lose'] = mne.read_evokeds(os.path.join(GApath, 'outcome_Lose_GA.fif'))[0]

filepath = os.path.join(evokedspath, 'outcome_Win.pkl')
with open(filepath, 'rb') as f:
    win = pickle.load(f)

filepath = os.path.join(evokedspath, 'outcome_Lose.pkl')
with open(filepath, 'rb') as f:
    lose = pickle.load(f)

for erp_type in erp_types:
    for sub in range(len(win)):

        lat_win[erp_type].append(extract_evoked_latency(win[sub], erp_type)) 
        lat_lose[erp_type].append(extract_evoked_latency(lose[sub], erp_type)) 

    if ttest_rel(np.array(lat_win[erp_type]),np.array(lat_lose[erp_type])).pvalue < significance_level_for_peak_differences:
        print('%s has different peak latencies between win and lose conditions'%erp_type)
        peak_latencies[erp_type+'_Win'] = round(extract_evoked_latency(GA['Win'], erp_type)[0], round_num)
        peak_latencies[erp_type+'_Lose'] = round(extract_evoked_latency(GA['Lose'], erp_type)[0], round_num)
    else: 
       peak_latencies[erp_type] = round(extract_evoked_latency(GA['overall'], erp_type)[0], round_num)

GA['sys'] = mne.read_evokeds(os.path.join(GApath, 'sysDiaFeedb_Systole_GA.fif'))[0]
GA['dia'] = mne.read_evokeds(os.path.join(GApath, 'sysDiaFeedb_Diastole_GA.fif'))[0]

filepath = os.path.join(evokedspath, 'sysDiaFeedb_Systole.pkl')
with open(filepath, 'rb') as f:
    sys = pickle.load(f)

filepath = os.path.join(evokedspath, 'sysDiaFeedb_Diastole.pkl')
with open(filepath, 'rb') as f:
    dia = pickle.load(f)

lat_sys = defaultdict(list)
lat_dia = defaultdict(list)

for erp_type in erp_types:
    for sub in range(len(sys)):

        lat_sys[erp_type].append(extract_evoked_latency(sys[sub], erp_type)) 
        lat_dia[erp_type].append(extract_evoked_latency(dia[sub], erp_type)) 

    if ttest_rel(np.array(lat_sys[erp_type]),np.array(lat_dia[erp_type])).pvalue < significance_level_for_peak_differences:
        print('%s has different peak latencies between systole and diastole conditions'%erp_type)
        peak_latencies[erp_type+'_Systole'] = round(extract_evoked_latency(GA['sys'], erp_type)[0], round_num)
        peak_latencies[erp_type+'_Diastole'] = round(extract_evoked_latency(GA['dia'], erp_type)[0], round_num)

if check_valence_phase_interactions:

    different_win_lose = defaultdict(list)
    different_sys_dia = defaultdict(list)

    # first, compare if win and lose are different within sys and dia separately
    for phase in ['Systole','Diastole']:
        GA['Win_' + phase] = mne.read_evokeds(os.path.join(GApath, 'sysDiaFeedb_' + phase + '_outcome_Win_GA.fif'))[0]
        GA['Lose_' + phase] = mne.read_evokeds(os.path.join(GApath, 'sysDiaFeedb_' + phase + '_outcome_Lose_GA.fif'))[0]
    
        filepath = os.path.join(evokedspath, 'sysDiaFeedb_' + phase + '_outcome_Win.pkl')
        with open(filepath, 'rb') as f:
            win = pickle.load(f)
        
        filepath = os.path.join(evokedspath, 'sysDiaFeedb_' + phase + '_outcome_Lose.pkl')
        with open(filepath, 'rb') as f:
            lose = pickle.load(f)
        
        for erp_type in erp_types:
            for sub in range(len(win)):
        
                lat_win[erp_type].append(extract_evoked_latency(win[sub], erp_type)) 
                lat_lose[erp_type].append(extract_evoked_latency(lose[sub], erp_type)) 

            different_win_lose[erp_type] = 0
            
            if ttest_rel(np.array(lat_win[erp_type]),np.array(lat_lose[erp_type])).pvalue < significance_level_for_peak_differences:
                print('%s has different peak latencies between win and lose conditions within %s' %(erp_type, phase))
                peak_latencies[erp_type+'_Win_'+ phase] = round(extract_evoked_latency(GA['Win_' + phase], erp_type)[0], round_num)
                peak_latencies[erp_type+'_Lose_'+ phase] = round(extract_evoked_latency(GA['Lose_' + phase], erp_type)[0], round_num)
                if abs(peak_latencies[erp_type+'_Win_'+ phase] - peak_latencies[erp_type+'_Lose_'+ phase]) > window_length/2:
                    different_win_lose[erp_type] = 1
                    print('%s peak latencies between win and lose conditions within %s are far enough' %(erp_type, phase))

    # then, can also compare if sys and dia are different within win and lose
    for valence in ['Win','Lose']:
        
        filepath = os.path.join(evokedspath, 'sysDiaFeedb_Systole_outcome_'+valence+'.pkl')
        with open(filepath, 'rb') as f:
            sys = pickle.load(f)
        
        filepath = os.path.join(evokedspath, 'sysDiaFeedb_Diastole_outcome_'+valence+'.pkl')
        with open(filepath, 'rb') as f:
            dia = pickle.load(f)
        
        lat_sys = defaultdict(list)
        lat_dia = defaultdict(list)
        
        for erp_type in erp_types:
            for sub in range(len(sys)):
        
                lat_sys[erp_type].append(extract_evoked_latency(sys[sub], erp_type)) 
                lat_dia[erp_type].append(extract_evoked_latency(dia[sub], erp_type)) 
            
            different_sys_dia[erp_type] = 0 
        
            if ttest_rel(np.array(lat_sys[erp_type]),np.array(lat_dia[erp_type])).pvalue < significance_level_for_peak_differences:
                print('%s has different peak latencies between systole and diastole within %s condition' %(erp_type, valence))
                peak_latencies[erp_type + '_' + valence + '_Systole'] = round(extract_evoked_latency(GA[valence + '_Systole'], erp_type)[0], round_num)
                peak_latencies[erp_type + '_' + valence + '_Diastole'] = round(extract_evoked_latency(GA[valence + '_Diastole'], erp_type)[0], round_num)
                if abs(peak_latencies[erp_type + '_' + valence + '_Systole'] - peak_latencies[erp_type + '_' + valence + '_Diastole']) > window_length/2:
                    different_sys_dia[erp_type] = 1
                    print('%s peak latencies between systole and diastole within %s condition are far enough' %(erp_type, valence))  

    for erp_type in erp_types:
        
        if different_sys_dia[erp_type] + different_win_lose[erp_type] == 0:
            print('for %s, take a peak latency based on average response across all trials' %erp_type)
            final_latencies[erp_type] = peak_latencies[erp_type]
            
        elif different_sys_dia[erp_type] + different_win_lose[erp_type] == 2:
            print('for %s, take peak latencies based on interactions of win-lose and systole-diastole conditions' %erp_type)
            final_latencies[erp_type + '_Win_Systole'] = peak_latencies[erp_type + '_Win_Systole']
            final_latencies[erp_type + '_Lose_Systole'] = peak_latencies[erp_type + '_Lose_Systole']
            final_latencies[erp_type + '_Win_Diastole'] = peak_latencies[erp_type + '_Win_Diastole']
            final_latencies[erp_type + '_Lose_Diastole'] = peak_latencies[erp_type + '_Lose_Diastole']
            
        elif different_win_lose[erp_type] == 1:
            print('for %s, take peak latencies based on win-lose conditions' %erp_type)
            final_latencies[erp_type + '_Win'] = peak_latencies[erp_type + '_Win']
            final_latencies[erp_type + '_Lose'] = peak_latencies[erp_type + '_Lose']
            
        else:
            print('for %s, take peak latencies based on systole-diastole conditions' %erp_type)
            final_latencies[erp_type + '_Systole'] = peak_latencies[erp_type + '_Systole']
            final_latencies[erp_type + '_Diastole'] = peak_latencies[erp_type + '_Diastole']

pd.DataFrame(peak_latencies, index=[0]).to_csv(os.path.join(derivativespath, 'all_significant_latencies.csv'), index=False)
if final_latencies:
    pd.DataFrame(final_latencies, index=[0]).to_csv(os.path.join(derivativespath, 'final_latencies.csv'), index=False)

evoked_pairs = [('Win', 'Lose'), ('sys', 'dia')]
evoked_pairs_interactions = [('Win_Systole', 'Lose_Systole'), ('Win_Diastole', 'Lose_Diastole')]

if save_vis:

    if not final_latencies:
        latencies = peak_latencies
    else:
        latencies = final_latencies
    
    for lat in latencies.keys():

        if lat.startswith('p3b'):
            channel = 'Pz'
        else:
            channel = 'FCz'

        ch_index = GA['overall'].ch_names.index(channel) 
        time_range_of_interest = [(latencies[lat] - window_length/2, latencies[lat] + window_length/2)]

        fig, axes = plt.subplots(3, 1, figsize=(10, 12))
        
        ax = axes[0]
        ax.plot(GA['overall'].times, GA['overall'].data[ch_index], label='overall', color='black')

        for start, end in time_range_of_interest:
            ax.axvspan(start, end, color='yellow', alpha=0.3)

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Amplitude (µV)')
        ax.legend()

        for i, (evoked_a, evoked_b) in enumerate(evoked_pairs):

            ax = axes[i+1]
        
            # Plot evoked_a
            ax.plot(GA[evoked_a].times, GA[evoked_a].data[ch_index], label=evoked_a, color='red')
            
            # Plot evoked_b
            ax.plot(GA[evoked_b].times, GA[evoked_b].data[ch_index], label=evoked_b, color='blue')
            
            # Highlight intervals of interest
            for start, end in time_range_of_interest:
                ax.axvspan(start, end, color='yellow', alpha=0.3)
            
            # Add labels and a title for each subplot
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Amplitude (µV)')
            ax.legend()

        plt.suptitle('Evoked Responses at %s Electrode, %s peak latency' %(channel,lat), fontsize=16)
        plt.tight_layout()
        fig.savefig(os.path.join(vispath, lat +'_evoked.jpg'), dpi=150, format='jpg')
        plt.close(fig)
        
        if final_latencies:

            fig, axes = plt.subplots(2, 1, figsize=(10, 12))
    
            for i, (evoked_a, evoked_b) in enumerate(evoked_pairs_interactions):
    
                ax = axes[i]
            
                # Plot evoked_a
                ax.plot(GA[evoked_a].times, GA[evoked_a].data[ch_index], label=evoked_a, color='red')
                
                # Plot evoked_b
                ax.plot(GA[evoked_b].times, GA[evoked_b].data[ch_index], label=evoked_b, color='blue')
                
                # Highlight intervals of interest
                for start, end in time_range_of_interest:
                    ax.axvspan(start, end, color='yellow', alpha=0.3)
                
                # Add labels and a title for each subplot
                ax.set_xlabel('Time (s)')
                ax.set_ylabel('Amplitude (µV)')
                ax.legend()
    
            plt.suptitle('Evoked Responses at %s Electrode, %s peak latency' %(channel,lat), fontsize=16)
            plt.tight_layout()
            fig.savefig(os.path.join(vispath, lat +'_evoked_interactions.jpg'), dpi=150, format='jpg')
            plt.close(fig)   

    time_ranges_of_interest = [lat for lat in latencies.values()]
    
    for key in GA.keys():
        
        fig = GA[key].plot_topomap(times=time_ranges_of_interest, sphere=(0.00, 0.02, 0.00, 0.09), average = window_length, show=False)
        fig.savefig(os.path.join(vispath, key +'_topo.jpg'), dpi=150, format='jpg')
        plt.close(fig)

    for (evoked_a, evoked_b) in evoked_pairs:

        dif = mne.combine_evoked([GA[evoked_a], GA[evoked_b]],[1,-1])
        fig = dif.plot_topomap(times=time_ranges_of_interest, sphere=(0.00, 0.02, 0.00, 0.09), average = window_length, show=False)
        fig.savefig(os.path.join(vispath, '%s_%s_topo.jpg'%(evoked_a, evoked_b)), dpi=150, format='jpg')
        plt.close(fig)

    if final_latencies:
        for (evoked_a, evoked_b) in evoked_pairs_interactions:

            dif = mne.combine_evoked([GA[evoked_a], GA[evoked_b]],[1,-1])
            fig = dif.plot_topomap(times=time_ranges_of_interest, sphere=(0.00, 0.02, 0.00, 0.09), average = window_length, show=False)
            fig.savefig(os.path.join(vispath, '%s_%s_topo.jpg'%(evoked_a, evoked_b)), dpi=150, format='jpg')
            plt.close(fig)
