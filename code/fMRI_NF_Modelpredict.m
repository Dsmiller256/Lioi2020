% Predicting fMRI-NF using EEG-NF

clear all;clc; 

% Loading Directories
cd ..
basedir = [pwd filesep];
savedir = [basedir, 'analysis' filesep 'eeg_analysis' filesep];
behdir = [basedir, 'analysis' filesep 'beh_analysis' filesep];
datadir = [basedir,'rawdata' filesep];
eegdir = [basedir, 'eeglab' filesep];

addpath(genpath(savedir));
addpath(genpath(basedir));
addpath(genpath(eegdir));
addpath(genpath(behdir));

%Load data
subjectlist={'sub-xp101','sub-xp102','sub-xp103','sub-xp104','sub-xp105','sub-xp106','sub-xp107','sub-xp108','sub-xp109','sub-xp110'};
fMRI_NF=[]; eeg_NF=[];
for sub=1:length(subjectlist)
    cd ([behdir subjectlist{sub} filesep 'NF_bold']);
    fNF=load(['d_' subjectlist{sub} '_task-fMRINF_NFbold_scores.mat']);
    eNF=load(['d_' subjectlist{sub} '_task-eegNF_NFbold_scores.mat']);
    
    fMRI_NF=[fMRI_NF;fNF.NF_bold.nf];
    eeg_NF=[eeg_NF;eNF.NF_bold.nf];
end

%Developing a model to predict fMRI activation

%The paper just says parameters without defining them. So I need to figure
%paraters to build the model!!