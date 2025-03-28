%% EEG preprocessing 

% Required dependencies:
% eeglab
% ANT EEProbe cnt reader plug-in
% HEPLAB plug-in for eeglab

% STEPS:
% Load data
% Find R-peaks and append to task events
% Save event data (behavior and Rpeaks)
% Low- and high-pass filter (0.1, 40) Hz
% Notch filter
% Remove bad channels, interpolate
% Re-reference to common average
% Save preprocessed files

close all;
clc; 
clear;

%% Define paths and settings

root = fileparts(pwd); 
data_path = [root '/data/'];
raw_path = [data_path 'raw/'];
preproc_path = [data_path 'preproc/'];
events_path = [data_path 'events/'];
misc_path = [data_path 'derivatives/'];

addpath(genpath([root '/codes/eeglab2024.1']));
addpath(genpath([root '/codes/ANTeepimport1.13']));
addpath(genpath([root '/codes/HEPLAB-master']));

eeglab;

% Set EEGLAB calculations to double precision
eeg_options
if option_single ~= 0 % This value should be zero
    option_single = 0;
end

rng(1234);

% If need to save additional files for further diagnostics 
save_files = true;

%% Find raw data files for automated preprocessing, and count subjects
rawDataFiles = dir([raw_path '*.cnt']);
subnum = length(rawDataFiles);

%% Loop through subjects

% Store number of bad channels
rejnum = zeros(subnum,1);
% Store bad channel names
rejchan = cell(subnum,1);

for isub=1:subnum
    
    loadName = rawDataFiles(isub).name;
    fileName = loadName(end-5:end-4); %extract only the subject number

    %% load
    % load file
    EEG = pop_loadeep_v4([raw_path loadName],'triggerfile', 'on');
    EEG.setname = fileName;
    % load channel locations
    EEG = pop_chanedit(EEG, 'load',[],'lookup','standard-10-5-cap385.elp');
    % drop the trigger channel; events are already saved in EEG.event
    % online reference - channel No 36, A1
    EEG = pop_select( EEG,'nochannel',{'TRIGGER', 'A1', 'A2'});
    % force the sampling rate to be 500 Hz (instead of automatically identified 499 Hz)
    EEG.srate = 500;

    if save_files
        pop_saveset(EEG,'filename',[fileName '_raw.set'],'filepath', preproc_path);
    end

    %% detect R-peaks
    % for some participants, invert ECG
    if ismember(fileName, {'02', '03', '04', '26'})
        EEG.data(64,:) = (EEG.data(64,:)*-1);
    end
    % find amplitudes and latencies of R-peaks
    ecg_raw = EEG.data(64,:); % pick ECG channel
    [qrs_amp_raw,qrs_i_raw]=heplab_pan_tompkin(double(ecg_raw), EEG.srate, 0);
    
    %% append the detected R-peaks into the EEG events struct 
    k = length(EEG.event); % number of existing events in the EEG structure
    evt = num2cell(qrs_i_raw);
    types = repmat({'Rpeak'},1, length(evt));
    [EEG.event(1,k+1:k+length(qrs_i_raw)).latency] = evt{:}; % assign latencies
    [EEG.event(1,k+1:k+length(qrs_i_raw)).type] = types{:}; % assign type (Rpeak)
    EEG = eeg_checkset(EEG,'eventconsistency'); % sorts all events by increasing latencies

    %% save all events (task and R-peaks)
    events_all = struct2table(EEG.event); %convert to table 
    writetable(events_all, [events_path, fileName '_allevents.xlsx'], "WriteMode","replacefile");

    %% filter, including ECG and EOG
    % 4th order butterworth filter
    [b,a] = butter(2,0.1/(EEG.srate/2),'high'); % high pass filter
    EEG.data = filtfilt(b,a,double(EEG.data)')'; 
    [c,d] = butter(2,40/(EEG.srate/2), 'low'); % low pass filter
    EEG.data = filtfilt(c,d,double(EEG.data)')';

    %% notch filter
    EEG = pop_eegfiltnew(EEG, 'locutoff',48,'hicutoff',52,'revfilt',1,'plotfreqz',0);

    if save_files
        pop_saveset(EEG,'filename',[fileName '_filter.set'],'filepath', preproc_path);
    end

    %% remove bad channels, interpolate    
    EEG_clean = pop_select(EEG,'nochannel',{'EOGH','EOGV','EKG1'});   
    all_chanlocs = EEG_clean.chanlocs;
    initial_ch_num = EEG_clean.nbchan;
  
    EEG_clean = clean_artifacts(EEG_clean, 'ChannelCriterion', 0.85, 'FlatlineCriterion', 5, 'BurstCriterion', 'off', 'WindowCriterion', 'off', 'Highpass', 'off');
    rejnum(isub) = initial_ch_num-EEG_clean.nbchan;   
    EEG.etc.clean_channel_num = EEG_clean.nbchan;
    
    if isfield(EEG_clean.etc, 'clean_channel_mask')
        EEG_clean = pop_interp(EEG_clean, all_chanlocs, 'spherical');    
        EEG.data(1:61,:) = EEG_clean.data(1:61,:);
        EEG.etc.clean_channel_mask = EEG_clean.etc.clean_channel_mask;
        ch_names = arrayfun(@(x) x.labels, EEG.chanlocs, 'UniformOutput', false);
        EEG.etc.interpolated_channels = ch_names(~EEG_clean.etc.clean_channel_mask);
        rejchan{isub} = EEG.etc.interpolated_channels; 
    end
    clear EEG_clean;

    %% common average 
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, [], 'exclude', [62, 63, 64]);
    EEG = pop_select(EEG,'nochannel',{'initialReference'});

    %% save preprocessed files
    pop_saveset(EEG,'filename',[fileName '_preproc.set'],'filepath', preproc_path);

end

writetable(array2table(rejnum), [misc_path, 'number_interpolated_channels.xlsx'], "WriteMode","replacefile")
writetable(cell2table(rejchan), [misc_path, 'interpolated_channels.xlsx'], "WriteMode","replacefile")

% quality check: number of removed channels should be less than 30
indices = find(rejnum > 30);

if ~isempty(indices)
    disp('Some subjects have too many bad channels')
    for i = 1:length(indices)
        disp(['File with too many bad channels: ' rawDataFiles(indices(i)).name])
    end
end