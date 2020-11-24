%% Setup

% Set dirs
cd ..;

basedir = [pwd filesep];
scriptdir = [basedir,'analysis' filesep 'eeg_analysis' filesep];
datadir = [basedir,'rawdata' filesep];
eegdir = [basedir, 'fieldtrip' filesep];

addpath(genpath(scriptdir));
addpath(genpath(basedir));
addpath ('D:\Lioi2020\Lioi2020\fieldtrip')

group = importdata([basedir,'sublist.txt']); % one subj per line
ft_defaults

% Bands to filter
all_bands = {'alpha'}; 
for subj = 1:length(group)
    cfg=[];
    cfg.dataset = ['D:\Lioi2020\Lioi2020\rawdata\' group{subj} '\eeg\' group{subj} '_task-eegNF_eeg.eeg'];
    trialdata=ft_preprocessing(cfg);
    
    cfg.viewmode="vertical"
    cfg=ft_databrowser(cfg,data)
    
    cfg.artfctdef.reject='complete';
    cleandata=ft_rejectartifact(cfg,data);
    
    cfg=[];
    cfg.channel = 'EEG';
    ic_data = ft_componentanalysis(cfg,cleandata)
    
    
end