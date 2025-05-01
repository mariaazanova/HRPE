import numpy as np
import pandas as pd
import mne
from mne import io 
import os
from support_functions import get_file_prefixes, find_next_same_cue_index, setup_logging, initialise_paths 
import logging

paths = initialise_paths()
root = paths['root']
eventspath = paths['eventspath']
cleanpath = paths['cleanpath']

logger = setup_logging(paths['logspath'], 'clearn06')

taskdatapath = os.path.join(root, 'CLEARN_longdata.xlsx')
taskdata = pd.read_excel(taskdatapath)
subjects = get_file_prefixes(cleanpath,'.fdt','.')

# load bad epochs/subjects from bad_ecg.csv (obtained at step 4) to mark them in df
# otherwise Rpeak event times should have already been changed in event files directly

logger.info('Started collecting behavioral and cardiac data')

df = pd.DataFrame(index=range(len(subjects)*300), columns=['subID', 
                                                           'block',
                                                           'miniBlock',
                                                           'half',
                                                           'trialCount', 
                                                           'trialInBlock',
                                                           'trialInMiniBlock',
                                                           'first25',
                                                           'next25',
                                                           'secondlast25',
                                                           'last25',
                                                           'learning_miniblock',
                                                           'first_next',
                                                           'first_reversal',
                                                           'first_last',
                                                           'reversal_next',
                                                           'reversal_last',
                                                           'next_last',
                                                           'RT',
                                                           'win',
                                                           'valence',
                                                           'outcome',
                                                           'choiceApproach',
                                                           'trueAversive',
                                                           'trueCorrect',
                                                           'trueCorrect_nextTrial',
                                                           'switchStay_nextTrialAction',
                                                           'switchNeg_nextTrialAction',
                                                           'stayPos_nextTrialAction',
                                                           'switchStay_sameCueAction',
                                                           'switchNeg_sameCueAction',
                                                           'stayPos_sameCueAction',
                                                           'switchStay_nextTrialStrategy',
                                                           'switchNeg_nextTrialStrategy',
                                                           'stayPos_nextTrialStrategy',
                                                           'percent_trueCorrect',
                                                           'percent_switchStay_nextTrialAction',
                                                           'percent_switchNeg_nextTrialAction',
                                                           'percent_stayPos_nextTrialAction',
                                                           'percent_switchStay_sameCueAction',
                                                           'percent_switchNeg_sameCueAction',
                                                           'percent_stayPos_sameCueAction',
                                                           'percent_switchStay_nextTrialStrategy',
                                                           'percent_switchNeg_nextTrialStrategy',
                                                           'percent_stayPos_nextTrialStrategy',
                                                           'RestAveIBI',
                                                           'RestRMSSD',
                                                           'RestSDNN',
                                                           'IBI1',
                                                           'IBI2',
                                                           'IBI3',
                                                           'IBI4',
                                                           'IBI5',
                                                           'IBI6',
                                                           'IBI7',
                                                           'IBI8',
                                                           'stimIBI1',
                                                           'stimIBI2',
                                                           'stimIBI3',
                                                           'stimIBI4',
                                                           'stimIBI5',
                                                           'stimIBI6',
                                                           'msDistToFeedb',
                                                           'msDistToStim',
                                                           'phaseFeedb',
                                                           'phaseStim',
                                                           'sysFeedb',
                                                           'sysStim',
                                                           'sysDiaFeedb',
                                                           'sysDiaStim',
                                                           'beforeFeedbSDNN',
                                                           'beforeFeedbRMSSD',
                                                           'afterFeedbSDNN',
                                                           'afterFeedbRMSSD',
                                                           'beforeStimSDNN',
                                                           'beforeStimRMSSD',
                                                           'afterStimSDNN',
                                                           'afterStimRMSSD',
                                                           'beforeFeedbAveIBI',
                                                           'afterFeedbAveIBI',
                                                           'beforeStimAveIBI',
                                                           'afterStimAveIBI',
                                                           'changeIBI2_1',
                                                           'changeIBI3_2',
                                                           'changeIBI4_3',
                                                           'changeIBI5_4',
                                                           'changeIBI6_5',
                                                           'changeIBI7_6',
                                                           'changeIBI8_7',
                                                           'changeStimIBI2_1',
                                                           'changeStimIBI3_2',
                                                           'changeStimIBI4_3',
                                                           'changeStimIBI5_4',
                                                           'changeStimIBI6_5'])

df['subID'] = np.repeat(subjects,300)
df['block'] = np.tile(np.repeat(list(range(1,4)),100),len(subjects))
df['miniBlock'] = np.tile(np.repeat(list(range(1,7)),50),len(subjects))
df['half'] = np.tile(np.repeat([1,2],50),3*len(subjects))
df['trialCount'] = np.tile(list(range(1,301)),len(subjects))
df['trialInBlock'] = np.tile(list(range(1,101)),3*len(subjects))
df['trialInMiniBlock'] = np.tile(list(range(1,51)),2*3*len(subjects))

taskdata.drop(taskdata[taskdata.subID == 17].index, inplace=True) 
df['RT'] = taskdata.RT.values
df['win'] = taskdata.correctResp.values
df['choiceApproach'] = taskdata.choiceApproach.values
df['trueAversive'] = taskdata.trueAversive.values
df['trueCorrect'] = taskdata.trueCorrect.values
del taskdata

df.loc[df.win==0,'valence'] = -1
df.loc[df.win==1,'valence'] = 1

df.loc[df.win==0,'outcome'] = 'Lose'
df.loc[df.win==1,'outcome'] = 'Win'

df['first25']=0
df.loc[df.trialInBlock<26,'first25']=1
df['next25']=0
df.loc[(df.trialInBlock>25)&(df.trialInBlock<51),'next25']=1
df['secondlast25']=0
df.loc[(df.trialInBlock>50)&(df.trialInBlock<76),'secondlast25']=1
df['last25']=0
df.loc[df.trialInBlock>75,'last25']=1

df['learning_miniblock'] = pd.cut(
    df['trialInBlock'],
    bins=[0, 25, 50, 75, 100],  # Define the edges for each category
    labels=['first', 'next', 'reversal', 'last'],
    right=True  # Include the right edge in each bin
)

df['first_next'] = df.apply(lambda row: 'first' if row['first25'] == 1 else ('next' if row['next25'] == 1 else None), axis=1)
df['first_reversal'] = df.apply(lambda row: 'first' if row['first25'] == 1 else ('reversal' if row['secondlast25'] == 1 else None), axis=1)
df['first_last'] = df.apply(lambda row: 'first' if row['first25'] == 1 else ('last' if row['last25'] == 1 else None), axis=1)
df['reversal_next'] = df.apply(lambda row: 'reversal' if row['secondlast25'] == 1 else ('next' if row['next25'] == 1 else None), axis=1)
df['reversal_last'] = df.apply(lambda row: 'reversal' if row['secondlast25'] == 1 else ('last' if row['last25'] == 1 else None), axis=1)
df['next_last'] = df.apply(lambda row: 'next' if row['next25'] == 1 else ('last' if row['last25'] == 1 else None), axis=1)

# Iterate through the df to populate the columns
for i in range(len(df)):
    next_index = i + 1 if i + 1 < len(df) else None
    
    # For next trial
    if next_index and df.loc[i, 'subID'] == df.loc[next_index, 'subID'] and df.loc[i, 'block'] == df.loc[next_index, 'block']:

        df.loc[i, 'trueCorrect_nextTrial'] = df.loc[next_index, 'trueCorrect']
        
        if df.loc[i, 'win'] == 1:
            df.loc[i, 'stayPos_nextTrialAction'] = int(df.loc[i, 'choiceApproach'] == df.loc[next_index, 'choiceApproach'])
            df.loc[i, 'stayPos_nextTrialStrategy'] = int(df.loc[i, 'trueCorrect'] == df.loc[next_index, 'trueCorrect'])
        elif df.loc[i, 'win'] == 0:
            df.loc[i, 'switchNeg_nextTrialAction'] = int(df.loc[i, 'choiceApproach'] != df.loc[next_index, 'choiceApproach'])
            df.loc[i, 'switchNeg_nextTrialStrategy'] = int(df.loc[i, 'trueCorrect'] != df.loc[next_index, 'trueCorrect'])
    
    # Finding next trial with same cue within the block
    next_same_cue_index = find_next_same_cue_index(df,i)
    if next_same_cue_index is not None:
        if df.loc[i, 'win'] == 1:
            df.loc[i, 'stayPos_sameCueAction'] = int(df.loc[i, 'choiceApproach'] == df.loc[next_same_cue_index, 'choiceApproach'])
        elif df.loc[i, 'win'] == 0:
            df.loc[i, 'switchNeg_sameCueAction'] = int(df.loc[i, 'choiceApproach'] != df.loc[next_same_cue_index, 'choiceApproach'])
            
df['switchStay_nextTrialAction'] = df['switchNeg_nextTrialAction'].combine_first(df['stayPos_nextTrialAction'])
df['switchStay_nextTrialStrategy'] = df['switchNeg_nextTrialStrategy'].combine_first(df['stayPos_nextTrialStrategy'])
df['switchStay_sameCueAction'] = df['switchNeg_sameCueAction'].combine_first(df['stayPos_sameCueAction'])

columns_to_compute_percentage = ['trueCorrect',
                                 'switchStay_nextTrialAction', 'switchNeg_nextTrialAction', 'stayPos_nextTrialAction',
                                 'switchStay_sameCueAction', 'switchNeg_sameCueAction', 'stayPos_sameCueAction',
                                 'switchStay_nextTrialStrategy', 'switchNeg_nextTrialStrategy', 'stayPos_nextTrialStrategy']

percent_columns = ['percent_' + column for column in columns_to_compute_percentage]

ibi_labels = ['IBI' + str(i) for i in range(1, 9)]
stim_ibi_labels = ['stimIBI' + str(i) for i in range(1, 7)]
change_labels = ['change' + ibi_labels[i] + '_' + str(i) for i in range(1,len(ibi_labels))]
change_stim_labels = ['change' + stim_ibi_labels[i].replace('s','S') + '_' + str(i) for i in range(1,len(stim_ibi_labels))]

for subject in subjects:

    logger.info(f"Processing subject {subject}")
    
    # load events
    subject_eventspath = os.path.join(eventspath, subject + '_allevents.xlsx')
    events = pd.read_excel(subject_eventspath)
    events.drop('duration', axis=1, inplace=True)
    events['latency']=events['latency']*2 # convert to ms, given 500 Hz sampling frequency
    
    # subject 01 is missing the first trigger for resting state
    if subject=='01':   
        events.loc[2,'type']='200'
        
    pick_sub = df['subID']==subject
    
    # for this subject, compute percent of some decision among all non-nan rows
    dfSub = df.loc[pick_sub, columns_to_compute_percentage]
    percent_ones = dfSub.apply(lambda col: (col == 1).sum() / col.notna().sum()).values
    df.loc[pick_sub, percent_columns] = percent_ones 
    
    # find latencies of events
    feedb_latencies = events.loc[((events.type=='51')|(events.type=='52')|(events.type=='53')|(events.type=='54')),'latency'].values
    stim_latencies = events.loc[((events.type=='21')|(events.type=='22')),'latency'].values
    Rpeak_latencies = events.loc[(events.type=='Rpeak'),'latency'].values
    Rest_latencies = events.loc[(events.type=='200'),'latency'].values

    Rest_Rpeaks = Rpeak_latencies[(Rpeak_latencies>Rest_latencies[0])&(Rpeak_latencies<Rest_latencies[1])]
    Rest_IBI = np.diff(Rest_Rpeaks)
    df.loc[pick_sub,'RestAveIBI'] = np.mean(Rest_IBI)
    df.loc[pick_sub,'RestSDNN'] = np.std(Rest_IBI, ddof=1)
    df.loc[pick_sub,'RestRMSSD'] = np.sqrt(np.mean(np.diff(Rest_IBI) ** 2))
    
    # pick 8 IBI surrounding feedback + find dist from feedback to previous Rpeak
    for feedb in range(len(feedb_latencies)):    
        lat_before_feedb = Rpeak_latencies[Rpeak_latencies<feedb_latencies[feedb]][-5:]
        lat_after_feedb = Rpeak_latencies[Rpeak_latencies>feedb_latencies[feedb]][:4]
        feedb_IBI = np.diff(np.concatenate([lat_before_feedb, lat_after_feedb]))
        df.loc[pick_sub&((df.trialCount-1)==feedb),ibi_labels] = feedb_IBI
        df.loc[pick_sub&((df.trialCount-1)==feedb),'msDistToFeedb'] = feedb_latencies[feedb]-lat_before_feedb[-1]        
    
    # 6 IBI surrounding stimulus presentation
    for stim in range(len(stim_latencies)):    
        lat_before_stim = Rpeak_latencies[Rpeak_latencies<stim_latencies[stim]][-4:]
        lat_after_stim = Rpeak_latencies[Rpeak_latencies>stim_latencies[stim]][:3]
        stim_IBI = np.diff(np.concatenate([lat_before_stim, lat_after_stim]))
        df.loc[pick_sub&((df.trialCount-1)==stim),stim_ibi_labels] = stim_IBI    
        df.loc[pick_sub&((df.trialCount-1)==stim),'msDistToStim'] = stim_latencies[stim]-lat_before_stim[-1]

df['phaseFeedb'] = (df.msDistToFeedb/df.IBI5)*360
df['phaseStim'] = (df.msDistToStim/df.stimIBI4)*360
df['sysFeedb'] = 0
df['sysStim'] = 0
df.loc[df.phaseFeedb<120,'sysFeedb'] = 1
df.loc[df.phaseStim<120,'sysStim'] = 1
df.loc[df.phaseFeedb<120,'sysDiaFeedb'] = 'Systole'
df.loc[df.phaseStim<120,'sysDiaStim'] = 'Systole'
df.loc[df.phaseFeedb>240,'sysDiaFeedb'] = 'Diastole'
df.loc[df.phaseStim>240,'sysDiaStim'] = 'Diastole'

before_feedb = ibi_labels[0:4]
before_stim = stim_ibi_labels[0:3]
after_feedb = ibi_labels[4:]
after_stim = stim_ibi_labels[3:]

df['beforeFeedbAveIBI'] = df[before_feedb].mean(axis=1).values
df['afterFeedbAveIBI'] = df[after_feedb].mean(axis=1).values
df['beforeStimAveIBI'] = df[before_stim].mean(axis=1).values
df['afterStimAveIBI'] = df[after_stim].mean(axis=1).values

df['beforeFeedbSDNN'] = np.std(df[before_feedb].values.astype('float64'), ddof=1,axis=1)
df['beforeFeedbRMSSD'] = np.sqrt(np.mean(np.diff(df[before_feedb].values.astype('float64'))**2,axis=1))
df['afterFeedbSDNN'] = np.std(df[after_feedb].values.astype('float64'), ddof=1,axis=1)
df['afterFeedbRMSSD'] = np.sqrt(np.mean(np.diff(df[after_feedb].values.astype('float64'))**2,axis=1))

df['beforeStimSDNN'] = np.std(df[before_stim].values.astype('float64'), ddof=1,axis=1)
df['beforeStimRMSSD'] = np.sqrt(np.mean(np.diff(df[before_stim].values.astype('float64'))**2,axis=1))
df['afterStimSDNN'] = np.std(df[after_stim].values.astype('float64'), ddof=1,axis=1)
df['afterStimRMSSD'] = np.sqrt(np.mean(np.diff(df[after_stim].values.astype('float64'))**2,axis=1))

df[change_labels] = np.diff(df[ibi_labels],axis=1)
df[change_stim_labels] = np.diff(df[stim_ibi_labels],axis=1)

df.to_csv(os.path.join(root, 'data.csv'), index = False)
logger.info('Saved single-trial data')

columns_of_interest = [col for col in df.columns if col.startswith("percent_") or col.startswith("Rest")]

unique_values_per_subject = []
for subID, group in df.groupby('subID'):
    unique_row = {'subID': subID}
    for col in columns_of_interest:
        unique_values = group[col].unique()
        if len(unique_values) > 1:
            logger.warning(f"Multiple unique values found for subject {subID} in column '{col}': {unique_values}")
        unique_row[col] = unique_values[0]  # Take the single unique value
    unique_values_per_subject.append(unique_row)

average_data = pd.DataFrame(unique_values_per_subject)
average_data.to_csv(os.path.join(root, "average_data.csv"), index=False)
logger.info('Saved subject-level average data')
