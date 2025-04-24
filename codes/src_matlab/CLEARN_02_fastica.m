%% Compute ICA 

% Required dependencies:
% eeglab
% fastICA

% STEPS:
% High-pass preprocessed continuous EEG at 1 Hz
% Crop from non-task activity
% Rank and run ICA (no component removal at this stage)
% Save ICA file

close all;
clc; 
clear;

%% Define paths and settings

root = fileparts(pwd); 
data_path = [root '/data/'];
ica_path = [data_path 'fastica/'];
preproc_path = [data_path 'preproc/'];

addpath(genpath([root '/codes/eeglab2024.1']));
addpath(genpath([root '/codes/FastICA_25']));

eeglab;

% Set EEGLAB calculations to double precision
eeg_options
if option_single ~= 0 % This value should be zero
    option_single =0;
end

rng(1234);

%% Find data files for automated processing, and count subjects
preprocFiles = dir([preproc_path '*_preproc.set']);
subnum = length(preprocFiles);

%% Loop through subjects

for isub=1:subnum
    
    loadName = preprocFiles(isub).name;
    tokens = regexp(loadName, '\d+', 'match'); % extract all number sequences
    fileName = tokens{1}; %extract only the subject number
    
    EEG = pop_loadset([fileName '_preproc.set'], preproc_path);   

    %% filter data for ICA with 1 Hz high-pass
    [b,a] = butter(2,1/(EEG.srate/2),'high'); 
    EEG.data = filtfilt(b,a,double(EEG.data)')'; 
    
    %% crop the recording
    
    % participant '01' does not have the first trigger for RS
    if ismember(fileName, {'01'})
        EEG.event(2).type = '200';
    end
    
    idx = find(contains({EEG.event.type}, {'200','101','102','103'})); % find beginnings and ends of 3 task blocks + resting-state
	breakpoints = [EEG.event(idx).latency];
	leave_points = [breakpoints(1) breakpoints(2); breakpoints(3) breakpoints(4); breakpoints(5) breakpoints(6); breakpoints(7) breakpoints(8)];

    EEG = pop_select(EEG, 'point', leave_points);

    %% rank 
    if isfield(EEG.etc, 'clean_channel_num')
        dataRank = EEG.etc.clean_channel_num; 
    else
        dataRank = rank(double(EEG.data'))-3; %-3 since we exclude ecg and eog electrodes in ica
    end

    %% run ICA
    EEG = pop_runica(EEG,'icatype','fastica','numOfIC', dataRank,'chanind',[1:61]);
    pop_saveset(EEG,'filename',[fileName '.set'],'filepath',ica_path); 

end
