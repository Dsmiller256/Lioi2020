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

%% fMRI Preprocessing (using SPM8)
spmdir = '/usr/local/MATLAB/spm8'; addpath(spmdir); spm fmri, close all

fs          = filesep; %makes code dynamic 

base_dir = '/mnt/datashare/Pr4.PF/Imaging/'; %base directory
data_dir = [base_dir 'Data']; %img data directory 
batch_dir = [base_dir 'Batch']; %batch files directory

%directories for different parts of the scan
dir_struct  = 'T1'; %anatomical scan
dir_fm      = 'Fieldmap';
sess_prfx   = 'R'; %EPI run: R1 R2 R3...
name_subj = {'
nRuns = 3;

%% Define what processing we want

fieldmap =1;%1;
biascorrect = 0; %skip bias correct unless you have additional sensors to correct movement
slicetime_correct =1;%1;
realign =1;%1;
coregister =1;%1;
segment_STRUCT =1;%1;
normalise =1;%1;
warpT1 =1;%1;
smooth =1;

%% MRI parameters: CHANGE ACCORDING TO YOUR MRI SEQUENCE!!
TR = 1.3; % TR in seconds
N_slices = 46; % number of slices
TA = TR-(TR/N_slices); % acquisition time for slice timing
RefSlice = 45; % reference slice for slice timing; middle slice
SliceOrder = [1:2:N_slices 2:2:N_slices];     % ascending slice order, 1:1:N_slices
FWHM  = 8; % smoothing kernel % was 6 in Marcos's script. will not use
% for now because doing MVPA
resolT1 = [1 1 1]; %resolution of structural (used during normalization)dfls
resolEPI = [3 3 3]; %resolution of functional (used during normalization)
% actual fx dims are 3x3x2.5
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corespondence subject/files
%for field maps

%-----------------------------------------------------------------------
% Slice time correction script
%-----------------------------------------------------------------------

for s0 = 1 : length(name_subj)

if fieldmap
disp(['Generating fieldmaps ']);
%for s0 = 1 : length(name_subj)

    disp(['Generating filedmaps for Subject : ', name_subj{s0}])

    fmDir = [data_dir fs name_subj{s0} fs dir_fm]; % go to the subject's fields dir
    fp   = spm_select('List', fmDir, ['^' name_subj{s0} '.*008_gre_field.*.img'] ); % Phase file (1 echo)
    fm   = spm_select('List', fmDir, ['^' name_subj{s0} '.*007.*01.img'] ); %magnitude file (2 echos)       

    phImage=cellstr([repmat([fmDir fs],size(fp,1),1) fp]);
    MagImage=cellstr([repmat([fmDir fs],size(fm,1),1) fm]);
    display(fp)
    display(fm)

    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase = phImage; %set phase img
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude = MagImage; %set mag img

    %one sub only got the "alex" field map so we have to apply a
    %different .m file
        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsfile = {'/usr/local/MATLAB/SPM/spm12/toolbox/FieldMap/pm_defaults.m'};


    for sess = 1:nRuns % if for all sessions seperate FMs then use -> 1:n_sess
        scanDir = [data_dir fs name_subj{s0} fs sess_prfx num2str(sess,'%01d')]; %second arugment in num2str is for number precision
        f = spm_select('List', scanDir, ['^' name_subj{s0} '.*_PS_UP30_X3_0006.img']);
        EPIfiles  = cellstr([repmat([scanDir fs],size(f,1),1) f]);
        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(sess).epi = EPIfiles;
        f = []; EPIfiles = [];
    end
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'session';
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 0;
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = ''; %AnatImage
    matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;

    %calculate vdm method was here?

    disp(['RUNNING field_maps ' name_subj{s0}]);
    spm_jobman('run',matlabbatch);
    clear matlabbatch
%end
end

if biascorrect
%for s0 = 1 : length(name_subj)
    disp(['Bias correction for Subject : ', name_subj{s0}]);
    for iR=1:nRuns
        epiDir = [data_dir fs name_subj{s0} fs 'run_' num2str(iR,'%01d')];
        spm_biascorrect(epiDir);
    end
%end
end


if slicetime_correct
%for s0 = 1 : length(name_subj)
    disp(['Slicetime correction for Subject : ', name_subj{s0}]);
    for sess = 1:nRuns
        scanDir = [data_dir fs name_subj{s0} fs sess_prfx num2str(sess,'%01d')];
        matlabbatch{1}.spm.temporal.st.scans{sess} = cfg_getfile('FPList', scanDir, ['^' name_subj{s0} '.*_PS_UP30_X3_.*.img']);
    end

    matlabbatch{1}.spm.temporal.st.nslices = N_slices;
    matlabbatch{1}.spm.temporal.st.tr = TR;
    matlabbatch{1}.spm.temporal.st.ta = TA;
    matlabbatch{1}.spm.temporal.st.so = SliceOrder;% 1:2:38 + 2:2:38;
    matlabbatch{1}.spm.temporal.st.refslice = RefSlice;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';

    disp(['%%%%%%%%% Starting slice time correction for subj: ', name_subj{s0},' %%%%%%%%%']);

    spm_jobman('run',matlabbatch);
    if ~exist([data_dir,'/',name_subj{s0},'/batchfiles_spm12/'])
        mkdir([data_dir,'/',name_subj{s0},'/batchfiles_spm12/'])
    end
    save([data_dir,'/',name_subj{s0},'/batchfiles_spm12/slice_time_correct.mat'],'matlabbatch')
    clear matlabbatch
%end
end



if realign
spm_figure('GetWin','Graphics'); % to make *.ps file
disp(['Realign: Realignment & UnWarp for MEAN and ALL images '])
%for s0 = 1 : length(name_subj)
    disp(['Realign job specification for Subject : ', name_subj{s0}]);

    % loop to define new session in job
    for sess = 1:nRuns
        %% find epi files
        scanDir = [data_dir fs name_subj{s0} fs sess_prfx num2str(sess,'%01d')];
        %f   = spm_select('List', scanDir, '^a.*_SSHB_up30_PA_X2_*.img');
        f   = spm_select('List', scanDir, '^a.*.img');
        files  = cellstr([repmat([scanDir '/'],size(f,1),1) f]);
        matlabbatch{1}.spm.spatial.realignunwarp.data(sess).scans = files;
        f = []; files = []; % clear temporary variables for next run
        fmDir = [data_dir fs name_subj{s0} fs dir_fm];
        fvdm5_files   = spm_select('List', fmDir, '^vdm5_.*img');
        fvdm5 = fvdm5_files(sess,:);
        vdm5Image = cellstr([repmat([fmDir fs],size(fvdm5,1),1) fvdm5]);
        matlabbatch{1}.spm.spatial.realignunwarp.data(sess).pmscan = vdm5Image;   % Note address the cell first with {}, then the structure with (), then put in cell... confusing!
    end

    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 2;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = {''};
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

    disp(['RUNNING realign and unwarp ' name_subj{s0}]);
    spm_jobman('run' , matlabbatch);
    save([data_dir,'/',name_subj{s0},'/batchfiles_spm12/realign.mat'],'matlabbatch')

%         newps=dir('*.ps');
%         isps=exist ('ps');
%         if isps~=7 %folder
%             mkdir('ps');
%         end
%         
%         for ips=1:numel(newps);
%             movefile(newps(ips).name, [pwd fs 'ps' fs name_subj{s0}, '_', newps(ips).name]);
%             delete(newps(ips).name);
%         end

    clear matlabbatch
end

if coregister
disp(['Coregister: Estimate structural (change postition file (.hdr) to match postition of T1w)'])
%for s0 = 1 : length(name_subj)
    disp(['Coregister job specification for Subject : ', name_subj{s0}])
    % first specify the structural as reference image
    stDir = [data_dir fs name_subj{s0} fs dir_struct];
    fT1  = spm_select('List', stDir, ['^' name_subj{s0} '.*_mprage_sag.*.img']);
    refImage = [stDir fs fT1];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref{:,1} = refImage;

    % then specify the first image as source image (can also use mean)
    sess=1;
    sourceDir = [data_dir fs name_subj{s0} fs sess_prfx num2str(sess,'%01d')];
    fmean  = spm_select('List', sourceDir, '^ua.*_PS_UP30_X3_0006.img$');
    sourceImage = [sourceDir fs fmean];
    matlabbatch{1}.spm.spatial.coreg.estimate.source{:,1} = sourceImage;

    % then specify other image, all EPIs
    otherImages=[]; f=[];
    for sess = 1:nRuns
        % find the epis
        scanDir = [data_dir fs name_subj{s0} fs sess_prfx num2str(sess, '%01d')];
        f  = spm_select('List', scanDir, '^ua*.*img');
        otherImages = [otherImages; cellstr([repmat([scanDir fs],size(f,1),1) f])];
    end

    for i = 1:size(otherImages,1)-1
        otherImages2{i,1} = otherImages{i+1,1};
    end

    matlabbatch{1}.spm.spatial.coreg.estimate.other = otherImages2;
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

    disp(['RUNNING coregister ' name_subj{s0}]);
    spm_jobman('run' , matlabbatch);
    save([data_dir,'/',name_subj{s0},'/batchfiles_spm12/coreg.mat'],'matlabbatch')
    clear matlabbatch
%end
end

if segment_STRUCT
disp(['Segmentation: Produce gray matter (native & modulated normalized) & white matter (native) from structural image']);
%for s0 = 1 : length(name_subj)
    TPM_file = [spmdir fs 'tpm' fs 'TPM.nii'];
    stDir = [data_dir fs name_subj{s0} fs dir_struct];

    fT1 = cfg_getfile('FPList', stDir, ['^' name_subj{s0} '.*_mprage_sag.*.img']);

    % Code from Mona
    preproc.channel.vols        = fT1;
    preproc.channel.biasreg     = 0.0001;   % bias regularisation cutoff
    preproc.channel.biasfwhm    = 60;       % bias FWHM
    preproc.channel.write       = [0 1];    % save bias corrected

    % tissue probability maps
    preproc.tissue(1).tpm       = {[TPM_file ',1']};
    preproc.tissue(1).ngaus     = 1;
    preproc.tissue(1).native    = [1 0];
    preproc.tissue(1).warped    = [0 0];
    preproc.tissue(2).tpm       = {[TPM_file ',2']};
    preproc.tissue(2).ngaus     = 1;
    preproc.tissue(2).native    = [1 0];
    preproc.tissue(2).warped    = [0 0];
    preproc.tissue(3).tpm       = {[TPM_file ',3']};
    preproc.tissue(3).ngaus     = 2;
    preproc.tissue(3).native    = [1 0];
    preproc.tissue(3).warped    = [0 0];
    preproc.tissue(4).tpm       = {[TPM_file ',4']};
    preproc.tissue(4).ngaus     = 3;
    preproc.tissue(4).native    = [1 0];
    preproc.tissue(4).warped    = [0 0];
    preproc.tissue(5).tpm       = {[TPM_file ',5']};
    preproc.tissue(5).ngaus     = 4;
    preproc.tissue(5).native    = [1 0];
    preproc.tissue(5).warped    = [0 0];
    preproc.tissue(6).tpm       = {[TPM_file ',6']};
    preproc.tissue(6).ngaus     = 2;
    preproc.tissue(6).native    = [0 0];
    preproc.tissue(6).warped    = [0 0];

    preproc.warp.mrf            = 1;
    preproc.warp.cleanup        = 0; % clean up any partitions: should get rid of unlikely structures (holes)
    preproc.warp.reg            = [0 0.001 0.5 0.05 0.2];
    preproc.warp.affreg         = 'mni';
    preproc.warp.fwhm           = 0;
    preproc.warp.samp           = 3;
    preproc.warp.write          = [0 1];

    matlabbatch{1}.spm.spatial.preproc = preproc;

    disp(['RUNNING segmemtation T1 images ' name_subj{s0}]);
    spm_jobman('run',matlabbatch);
    clear matlabbatch
%end
end

if normalise
disp(['Normalizing ... ']);
%for s0 = 1 : length(name_subj)
    wdir = [data_dir fs name_subj{s0} fs dir_struct];

    % Loop over sessions for epi's
    conCat_files=[];
    for sess = 1:nRuns
        scanDir = [data_dir fs name_subj{s0} fs sess_prfx num2str(sess,'%01d')];
        % select slice time corrected images
        f   = spm_select('List', scanDir, '^ua.*_PS_UP30_X3_.*.img');
        files  = cellstr([repmat([scanDir fs],size(f,1),1) f]);
        conCat_files = [conCat_files; files];       % concatenate all files over runs
        f = []; files = [];
    end


    fy  = spm_select('List', wdir, ['^y_' name_subj{s0} '.*.nii']); %write nomalize data
    yfile  = {[wdir fs fy]};
    matlabbatch{1}.spm.spatial.normalise.write.subj.vol = conCat_files(1,:);
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = yfile; % parameter file from segment
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = conCat_files;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -70; 78 76 85]; % Try changing to 2 2 2
    matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = resolEPI; % Try changing to 3 3 3
    matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 7; % changed from 4
    matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';
    disp(['RUNNING normalization' name_subj{s0}]);
    spm_jobman('run',matlabbatch);
    save([data_dir,'/',name_subj{s0},'/batchfiles_spm12/normalise.mat'],'matlabbatch')
    clear matlabbatch
%end
end

if warpT1
disp(['Normalizing T1s... ']);
%for s0 = 1 : length(name_subj)
    wdir = [data_dir fs name_subj{s0} fs dir_struct];

    %%% for the structurals
    fy  = spm_select('List', wdir, '^y_*.*');
    yfile  = {[wdir fs fy]};
    f = []; files = []; conCat_files = [];
    f  = spm_select('List', wdir, '^m.*_mprage_sag.*.nii');        
    files  = cellstr([repmat([wdir fs],size(f,1),1) f]); conCat_files = [conCat_files; files];

    matlabbatch{1}.spm.spatial.normalise.write.subj.def = yfile; % parameter file from segment
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = conCat_files;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.preserve = 0;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -70; 78 76 85];
    % use resolution 1 1 1 for structurals (i.e. T1w,MTw), but use 2 2 2 for EPI (functionals)
    matlabbatch{1}.spm.spatial.normalise.write.roptions.vox = resolT1;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.normalise.write.roptions.prefix = 'w';

    disp(['RUNNING normalization for T1s ' name_subj{s0}]);
    spm_jobman('run',matlabbatch);
    clear matlabbatch
%end
end

if smooth
disp(['Smoothing EPIs... ']);
epiimgs=[];
    for sess = 1:nRuns    
        if sess==1
            epiimgs2add = 0;
        else
            epiimgs2add = epiimgs2add+length(epiimgs);
        end            
        scanDir = [data_dir fs name_subj{s0} fs sess_prfx num2str(sess,'%2d')];
        epiimgs = spm_select('List',scanDir,'^wua.*.img'); %changed from wuabf %take off ed when done

        for epi = 1:length(epiimgs)
            spatial.smooth.data{epiimgs2add+epi,1} = [scanDir fs epiimgs(epi,:)];
            %spatial.smooth.data{1,1} = [scanDir fs epiimgs(epi,:)];
        end
    end

    spatial.smooth.fwhm     = [8 8 8];
    spatial.smooth.dtype    = 0;
    spatial.smooth.im       = 0;
    spatial.smooth.prefix   = 's';

    matlabbatch{1}.spm.spatial = spatial;
    disp(['%%%%%%%%% Starting smoothing for subj: ', name_subj{s0},' %%%%%%%%%']);
    spm_jobman('run',matlabbatch);
    clear matlabbatch;
end
end


% To include once I've sorted out how to do the physio analysis with Eric
% if physio
%     myPhysio_regressors_Hversion;
% end

        
        
        
 