%% Reject components

% Required dependencies:
% eeglab
% SASICA plugin

% STEPS:
% Load the initial 0.1 Hz dataset and ICA results on 1 hz dataset
% Run SASICA to automatically find components for rejection 
% Run ICLabel to identify potential components for rejection
% Combine results of SASICA and ICLabel
% Apply ICA results to 0.1 Hz high-passed dataset
% Save clean continuous data

close all;
clc;
clear;

%% Define paths and settings

root = fileparts(pwd); 
data_path = [root '/data/'];
ica_path = [data_path 'fastica/'];
preproc_path = [data_path 'preproc/'];
clean_path = [data_path 'clean/'];
vis_path = [root '/vis/rejectedComp/'];
misc_path = [data_path 'derivatives/'];

addpath(genpath([root '/codes/eeglab2024.1']));
addpath(genpath([root '/codes/SASICA-master']));

eeglab;

% Set EEGLAB calculations to double precision
eeg_options
if option_single ~= 0 % This value should be zero
    option_single =0;
end

rng(1234);

% If need to save plots of rejected components
save_plots = true;

%% Find ICA data files for automated processing, and count subjects
icaFiles = dir([ica_path '*.set']);
subnum = length(icaFiles);

%% Loop to identify and remove components
numComp = zeros(subnum,1); % track number of components
numRejComp = zeros(subnum,1); % track number of rejected components
numEcgComp = zeros(subnum,1); % track number of rejected components related to cardiac artifacts
RejComp = cell(subnum,1); % track indices of rejected components
EcgComp = cell(subnum,1); % track indices of rejected components related to cardiac artifacts

for isub=1:subnum

    loadName = icaFiles(isub).name;
    tokens = regexp(loadName, '\d+', 'match'); % extract all number sequences
    fileName = tokens{1}; %extract only the subject number

    EEG_ICA = pop_loadset([fileName '.set'], ica_path);  % ICA result on 1 Hz data
    EEG = pop_loadset([fileName '_preproc.set'], preproc_path); % 0.1 Hz preprocessed data

    %% use ICLabel
    EEG_ICA = iclabel(EEG_ICA); 

    %% use SASICA 
    close all
    cfg.opts.noplot = 1;
    cfg.opts.nocompute = 0;
    cfg.autocorr.enable = true;
    cfg.autocorr.dropautocorr = 'auto 2'; % threshold is determined automatically based on SD
    cfg.autocorr.autocorrint = 20;% will compute autocorrelation with this many milliseconds lag
    cfg.focalcomp.enable = true;
    cfg.focalcomp.focalICAout = 'auto 2'; % using default, which is 2SD - threshold is automatically determined
    cfg.trialfoc.enable = false;
    cfg.trialfoc.focaltrialout = 'auto 2'; % using default which is 2SD
    cfg.resvar.enable = false;
    cfg.resvar.thresh = 40 ;% %residual variance allowed, default is 15%
    cfg.SNR.enable = false;

    cfg.EOGcorr.enable = true;
    cfg.EOGcorr.corthreshV ='auto 4';% threshold correlation with vertical EOG, default would be 4 SD from average correlation: 'auto 4'
    cfg.EOGcorr.Veogchannames = 'EOGV';% vertical channel(s)
    cfg.EOGcorr.corthreshH ='auto 4';% threshold correlation with horizontal EOG
    cfg.EOGcorr.Heogchannames = 'EOGH'; %
    cfg.chancorr.enable = true;
    cfg.chancorr.channames ='EKG1';
    cfg.chancorr.corthresh='auto 2';

    cfg.FASTER.enable = false; % default thresh is 3 SDs from average for each measure
    cfg.FASTER.blinkchans = [];
    cfg.ADJUST.enable = false; % combines different features for detection (SAD,SVD,TK)
    cfg.MARA.enable = false;
    [EEG_ICA, cfg] = eeg_SASICA(EEG_ICA,cfg);

    %% Identify components to reject based on SASICA and ICLabel results
    
    % ICLabel results
    ICLabel_ECG_prob = EEG_ICA.etc.ic_classification.ICLabel.classifications(:,4); % prob. of ECG artifacts
    ICLabel_brain_prob = EEG_ICA.etc.ic_classification.ICLabel.classifications(:,1); % prob. of brain activity
    ICLabel_muscle_prob = EEG_ICA.etc.ic_classification.ICLabel.classifications(:,2); % prob. of muscle noise
    ICLabel_eye_prob = EEG_ICA.etc.ic_classification.ICLabel.classifications(:,3); % prob. of EOG artifacts
    ICLabel_other_max_prob = max(EEG_ICA.etc.ic_classification.ICLabel.classifications(:,5:7).').'; % prob. of other artifacts
    % Pick components that are likely artifactual given that they don't contain high probability of brain activity 
    brain_prob_thresh = 0.05;
    ECG_prob_thresh = 0.9;
    muscle_prob_thresh = 0.9;
    eye_prob_thresh = 0.9;
    other_prob_thresh = 0.9;
    ICLabel_ECG_comp = find((ICLabel_ECG_prob > ECG_prob_thresh)&(ICLabel_brain_prob < brain_prob_thresh)).';
    ICLabel_muscle_comp = find((ICLabel_muscle_prob > muscle_prob_thresh)&(ICLabel_brain_prob < brain_prob_thresh)).';
    ICLabel_eye_comp = find((ICLabel_eye_prob > eye_prob_thresh)&(ICLabel_brain_prob < brain_prob_thresh)).';
    ICLabel_other_comp = find((ICLabel_other_max_prob > other_prob_thresh)&(ICLabel_brain_prob < brain_prob_thresh)).';
    all_ICLabel_comp = unique(sort([ICLabel_ECG_comp ICLabel_muscle_comp ICLabel_eye_comp ICLabel_other_comp]));
    
    % SASICA components (+check that ICLabel's brain probability for these components is not high)
    all_SASICA_comp = find((EEG_ICA.reject.gcompreject == 1)&((ICLabel_brain_prob < brain_prob_thresh)')); 
    ECG_corr = EEG_ICA.reject.SASICA.icachancorrchans; % correlations of all components with ECG
    ECG_thresh = EEG_ICA.reject.SASICA.icathreshchancorrchans; % threshold for rejection of ECG components
    SASICA_ECG_comp = find((ECG_corr > ECG_thresh)&((ICLabel_brain_prob < brain_prob_thresh)')); 

    % combine SASICA and ICLabel results
    comprej = unique(sort([all_SASICA_comp all_ICLabel_comp])); % all components to be rejected
    ecgrej = unique(sort([SASICA_ECG_comp ICLabel_ECG_comp])); % ECG components to be rejected

    % track rejected componentns
    numComp(isub) = length(EEG_ICA.reject.gcompreject); % number of components 
    numRejComp(isub) = length(comprej); % number of rejected components
    numEcgComp(isub) = length(ecgrej); % number of rejected components related to cardiac artifacts
    RejComp{isub} = comprej; % indices of rejected components
    EcgComp{isub} = ecgrej; % indices of rejected components related to cardiac artifacts
    
    % plot components for rejection, if needed
    if save_plots       
        for icomp=1:length(comprej)
            figure; pop_prop_extended(EEG_ICA, 0, comprej(icomp));
            saveas(gcf,strcat(vis_path,fileName,'_rejected_component_',string(icomp),'.png'));
            close all;
        end 
    end

    %% transfer ICA result to 0.1 Hz dataset
    
    EEG.comprej = comprej;
    EEG.ecgrej = ecgrej;

    EEG.icawinv = EEG_ICA.icawinv;
    EEG.icasphere = EEG_ICA.icasphere;
    EEG.icaweights = EEG_ICA.icaweights;
    EEG.icachansind = EEG_ICA.icachansind;
    eeg_checkset(EEG);
        
    clear EEG_ICA;  

    %% reject chosen ICA components
    EEG = pop_subcomp(EEG, [EEG.comprej], 0); 

    %% save cleaned data
    pop_saveset(EEG,[fileName '.set'], clean_path);

end

%% save info about components
writetable(array2table(numComp), [misc_path, 'numComp.xlsx'], "WriteMode","replacefile")
writetable(array2table(numRejComp), [misc_path, 'numRejComp.xlsx'], "WriteMode","replacefile")
writetable(array2table(numEcgComp), [misc_path, 'numEcgComp.xlsx'], "WriteMode","replacefile")
writetable(cell2table(RejComp), [misc_path, 'RejComp.xlsx'], "WriteMode","replacefile")
writetable(cell2table(EcgComp), [misc_path, 'EcgComp.xlsx'], "WriteMode","replacefile")

% quality check: number of removed heart components on average is higher
% than 0.5
if mean(numEcgComp) < 0.5
    disp('Number of removed cardiac components is too low')
end

% quality check: for all subjects, at least one component was removed
indices_Comp = find(numComp < 1);
if ~isempty(indices_Comp)
    disp('No components were removed for some subjects')
    for i = 1:length(indices_Comp)
        disp(['Subjects with no components: ' icaFiles(indices_Comp(i)).name])
    end
end
