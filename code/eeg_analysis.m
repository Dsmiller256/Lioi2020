function eeg_analysis(step)
% HPF 0.5
% LPF 30
%
% Possible preprocessing steps:
%                               'preproc'

% Run with EEGLAB v. 14.0.0

%% Setup

% Set dirs
cd ..;

basedir = [pwd filesep];
scriptdir = [basedir,'analysis' filesep 'eeg_analysis' filesep];
datadir = [basedir,'rawdata' filesep];
eegdir = [basedir, 'eeglab' filesep];

addpath(genpath(scriptdir));
addpath(genpath(basedir));
addpath(genpath(eegdir));


group = importdata([basedir,'sublist.txt']); % one subj per line


% Bands to filter
all_bands = {'alpha','theta','beta','gamma'}; 

% initialize error log
errorlog = {}; ctr=1;

for subj = 1:length(group)
    try
        % Start EEGLAB
        cd (eegdir);
        [ALLEEG,EEG,CURRENTSET,ALLCOM] = eeglab;
        
        suffix = '_data';
        
        %% set up subj info
        subid = group{subj};
        fprintf(['\n\n',subid,'\n']);
        
        filepath = [datadir, subid, filesep 'eeg'];
        
        % cd into subj folder
        cd(filepath);
            runs = [1];
        end
        
        n_runs = size(runs,2);
        
        %% PREPROC STEP
        
        if strcmp(step, 'preproc')
            
            % Create vars for later channel rejection
            reject_indices = {};
            all_reject = [];
            
            for i = 1:1
                r = 1;
                 run_file = [datadir,subid,filesep, 'eeg' filesep subid '_task-eegNF_eeg.eeg'];
                 EEG = pop_readegi(run_file, [],[],'auto');
 
                curr_suffix = suffix;
                
                %% Load channel locations
                EEG = pop_chanedit(EEG,'load',{[basedir, 'analysis' filesep 'eeg_analysis' filesep 'GSN-HydroCel-128.sfp'] 'filetype' 'autodetect'});
                
                %% Remove noise at start of run
                EEG = auto_reject(EEG); % from impedance measurement
                
                %% Decimate from 1000Hz to 100Hz
                EEG = pop_resample(EEG, 200);
                
                %% Bandpass filter
                fprintf(['\nFiltering subject ',subid,', run ',num2str(r),'\n']);
                
                % Hamming windowed since FIR filter
                EEG = pop_eegfiltnew(EEG,[],50); % low pass filter
                
                suffix = [suffix, 'f'];
                
                % Save data
                [ALLEEG,EEG,CURRENTSET] = pop_newset(ALLEEG,EEG,1,'setname', ...
                    ['s',num2str(r),suffix],'savenew',['s',num2str(r), ...
                    suffix,'.set'],'overwrite','on','gui','off');
                
                %% Find bad channels to remove
                %EEG = pop_loadset('filename',['s',num2str(r),suffix,'.set']);
                
                [~,indelec] = pop_rejchanspec( EEG, 'stdthresh', [-3, 3]);
               
                reject_indices{r,1} = indelec;
                all_reject = [all_reject, reject_indices{r,1}];
                
                suffix = curr_suffix; % Reset suffix for next run
            end

            suffix = [suffix,'bf'];
            
            %% Remove bad channels from all runs
            for i = 1:1
                r = num2str(runs(i));
                
                curr_suffix = suffix;
                file = ['s',r,suffix,'.set'];
                EEG = pop_loadset('filename',file,'filepath',filepath);
                EEG = pop_select(EEG, 'nochannel', all_reject);
                suffix = [suffix,'b'];
                
                % Save out indices -- reject all for each run (union)
                save('rejected_channels', 'reject_indices');
                
                %% Save out data
                [ALLEEG,EEG,CURRENTSET] = pop_newset(ALLEEG,EEG,1,'setname', ...
                    ['s',num2str(r),suffix],'savenew',['s',r, ...
                    suffix,'.set'],'overwrite','on','gui','off');
                suffix = curr_suffix; % Reset suffix for next run
            end
        end
 