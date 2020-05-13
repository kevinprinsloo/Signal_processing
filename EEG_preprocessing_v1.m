
% Summary:
% Script to Analyse autotune MAMTA FS1 Data

% Status:
% Under Development

% Notes:
% n/a

% Author(s):
% Kevin Prinsloo

% Editor(s):
%

%% Prepare Workspace
%clearvars
%close all
%clc

%% Prepare variables for Cluster
% Manually Initialise Variables

% Initialise Path Variables
addpath '/scratch/kprinslo/Chimera_Study/Toolboxes/eeglab_current/eeglab14_1_2b'
eeglab
close all

addpath '/scratch/kprinslo/RI_Study/Toolboxes/fieldtrip_new/';
addpath '/scratch/kprinslo/RI_Study/Resources_Misc/';
addpath '/scratch/kprinslo/RI_Study/Scripts/';

% Setup Fiedltrip
ft_defaults

% Define study folder
study_name = 'IAPS_Pos_Neu_35';

% Initialise Subject Variables
listing = dir([study_path,'/',study_name,'/','BDF_data','/']);
subejct_listings = {listing.name};
subejct_listings(cellfun('length',subejct_listings)<3) = [];
subjects_orig = subejct_listings;
subjects_number = numel(subjects_orig);
subjects = subjects_orig;

% Initialise Condition Variables
triggers = [75,76,85,86,100,125,10,11,20,21,30,31,40,41];
%% Resp
conditions{1} = 'resp';
%% GO
conditions{75} = 'N_go';
conditions{76} = 'P_go';
conditions{10} = 'P_go';
conditions{20} = 'N_go';
conditions{30} = 'Cocaine_go';
conditions{40} = 'xy_go';
conditions{100} = 'P_go_100';
% Concat
triggers_Go = [75 76 10 20 30 40 100];

%% NoGo
conditions{85} = 'N_NoGo';
conditions{86} = 'P_NoGo';
conditions{11} = 'P_NoGo';
conditions{21} = 'N_NoGo';
conditions{31} = 'Cocaine_NoGo';
conditions{41} = 'xy_NoGo';
conditions{125} = 'P_NoGo_125';
% Concat
triggers_NoGo = [85 86 11 21 31 41 125];

%FILTER
FstopH = .0025;
FpassH = .01;
AstopH = 60;
FpassL = 40;
FstopL = 50;
AstopL = 40;
Apass = 1;

fs = 512;
Fs = 512;
% Generate high/low-pass filters
h = fdesign.highpass(FstopH,FpassH,AstopH,Apass,fs);
hpf = design(h,'cheby2','MatchExactly','stopband'); clear h
h = fdesign.lowpass(FpassL,FstopL,Apass,AstopL,fs);
lpf = design(h,'cheby2','MatchExactly','stopband'); clear h

% Initialise Stimulus Variables
stimulus_sampling_rate_original_Hz = 512;
stimulus_sampling_rate_downsampled = 256;

prestim = 0.1;
poststim = 0.8;

%--------------------------
%% Perform TRF analyses
%--------------------------

% Cluster parallel definition
subject_idx_cluster = str2double(getenv('SLURM_ARRAY_TASK_ID'));  % this give's you back the job parameter from slurm (#SBATCH --array=1-16)
disp(subject_idx_cluster)

% Look over each subject
for subject_idx = subject_idx_cluster
    % subject_idx = 1;
    
    % List BSF files
    clear data_pre
    % study_path_bdf_load = 'Z:\Neurotypical_Response_Inhibition\';
    listing = dir(fullfile([study_path,'/',study_name,'/','BDF_data','/',num2str(subjects_orig{subject_idx}),'/'], '*.bdf'));
    
    % Loop across BDF blocks
    data_ReRef = cell(length(listing),1);
    data_CSD_spline = cell(length(listing),1);
    data_rm_ds = cell(length(listing),1);
    for block_idx = 1:length(listing)
        
        % Define data file
        fullname = ([study_path,'/',study_name,'/','BDF_data','/',num2str(subjects_orig{subject_idx}),'/',listing(block_idx).name]); % EEG FILE
        
        % Access events and index
        cfg = [];
        cfg.dataset = fullname; % EEG data.
        event = ft_read_event(cfg.dataset);
        hdr   = ft_read_header(cfg.dataset);
        eventS = ft_filter_event(event, 'type','STATUS'); % find all triggers
        
        %% Epoching & Filtering
        % Epoch one continous dataset & apply filters
        cfg = [];
        cfg.headerfile = fullname;
        cfg.datafile = fullname;
        cfg.trialfun = 'ft_trialfun_general';
        cfg.channel = 'EEG';
        cfg.detrend = 'no';
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 45;
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.1;
        cfg.hpfilttype = 'firws';
        cfg.hpfiltdir = 'onepass-zerophase';
        cfg.hpfiltwintype = 'hamming';
        cfg.trialdef.triallength = Inf;
        cfg.trialdef.ntrials = 1;
        cfg = ft_definetrial(cfg);
        all_data = ft_preprocessing(cfg);
        
        %     % Filter EEG
        %     data_hold = all_data.trial{1}';
        %     data_hold = filtfilthd(hpf,data_hold);
        %     data_hold = filtfilthd(lpf,data_hold);
        %     all_data.trial{1} = data_hold';
        
        %% Define triggers for trialinfo
        holder = zeros(1,length(eventS));
        for k = 1:length(eventS)
            holder(k) =  eventS(k).value;
        end
        holder = holder';
        fprintf('trigger holder found\n');
        
        % Redefine required triggers
        trigRange = unique(holder);
        trigsNum = length(trigRange);
        if ~isequal(trigRange,triggers)
            cnt = setdiff(trigRange, triggers);
            [tf,idx] = ismember(trigRange, triggers);
            c_r = find(tf == 0);
            trigRange(c_r) =[];
            trigsNum = length(trigRange);
        end
        trigRange_resp = ([trigRange' 1]);
        
        % search for "trigger" events
        clear val t_onset value sample value_hld c
        c=1; % counter
        for k=1:length(eventS)
            if sum(eventS(k).value==trigRange_resp) && (eventS(k).sample>1)
                val(c) = eventS(k).value;
                t_onset(c) = eventS(k).sample;
                c=c+1;
            end
        end
        value  = val';
        sample = t_onset';
        value_hld = value;
        
        % determine the number of samples before and after the trigger
        pretrig = -prestim * hdr.Fs;
        posttrig = poststim * hdr.Fs;
        rsp_indx = find(ismember(value,1));
        stim_indx = find(ismember(value,trigRange_resp));
        trl = []; cnt2 = 1; not_inc = [];
        for istim = 2:length(stim_indx)
            if sample(istim) < abs(pretrig)
                fprintf('not inclucinf trial')
                value_hld(cnt2) = [];
                not_inc(cnt2) = istim;
                cnt2=cnt2+1;
            else
                if ismember(value(stim_indx(istim)), [trigRange_resp]) && istim < length(stim_indx) && (eventS(istim).sample>1)
                    trl_rsp_indx = find(rsp_indx > stim_indx(istim),1,'first');
                    rsp = value(rsp_indx(trl_rsp_indx));
                    RT  = (sample(rsp_indx(trl_rsp_indx)) - sample(stim_indx(istim))) / hdr.Fs;
                else
                    rsp = 0;
                    RT  = 0;
                end
                newtrl = [ round(sample(stim_indx(istim))+pretrig) floor(sample(stim_indx(istim))+posttrig)-1 round(pretrig) value(stim_indx(istim)) rsp RT];
                trl = [ trl; newtrl ];
            end
        end
        trl(end,:) = [];
        
        %% Epoch data
        clear dataX_epoch
        cfg = [];
        cfg.headerfile = fullname;  % name of EEG file
        cfg.dataset = fullname; % name of eeg file
        cfg.channel = 'EEG';
        cfg.continuous = 'yes';
        cfg.trl = trl;
        %cfg.trialfun = 'ft_trialfun_general';
        %cfg.trialdef.eventvalue = trigRange;
        %cfg.trialdef.prestim   = 0.3; % in seconds
        %cfg.trialdef.poststim  = 0.9; % in seconds
        %cfg.trialdef.eventtype = 'STATUS'; % get a list of the available types
        cfg = ft_definetrial(cfg);
        dataX_epoch = ft_redefinetrial(cfg,all_data); %redefines the filtered data
        clear all_data
        
        cfg = [];
        cfg.demean = 'yes';
        cfg.continuous = 'yes';
        dataX_epoch = ft_preprocessing(cfg,dataX_epoch);
        
        %% Define trialinfo
        trl_holder = value_hld(2:end-1);
        trl_holder_tm = zeros(length(trl_holder),2);
        for k = 3:length(trl_holder)-4            
            % Define range
            trgC = trl_holder(k-2:k+3);
            %% >>> Correct Rejection            
            if double(any(trgC(1) == triggers_Go)) && trgC(2) == 1 && double(any(trgC(3) == triggers_NoGo)) && trgC(4) ~= 1; trl_holder_tm(k,1) = 2; end
            %% >>> False Alarm
            if double(any(trgC(1) == triggers_Go)) && trgC(2) == 1 && double(any(trgC(3) == triggers_NoGo)) && trgC(4) == 1; trl_holder_tm(k,1) = -1; end
            %% >>> Hit
            if double(any(trgC(1) == triggers_Go)) && trgC(2) == 1 && double(any(trgC(3) == triggers_Go)) && trgC(4) == 1; trl_holder_tm(k,1) = 1; end
            %% >>> False Alarm Resp
            if double(any(trgC(1) == triggers_Go)) && trgC(2) == 1 && double(any(trgC(3) == triggers_NoGo)) && trgC(4) == 1; trl_holder_tm(k+1,2) = -1; end
        end
                   
        trl(:,5) = trl_holder_tm(:,1);
        trl(:,7) = trl_holder_tm(:,2);
        dataX_epoch.trialinfo = trl;
        
        %-----------------------
        %% Remove bad channels
        %-----------------------
        channels_number_cephalic = size(dataX_epoch.label,1);
        channel_locations = readlocs([study_path,'/','Resources_Misc','/','BioSemi','_',num2str(channels_number_cephalic),'_','AB','.sfp'],'filetype','sfp');
        
        clear dat_temp
        for j = 1:length(dataX_epoch.trial)
            dat_temp(:,:,j) = dataX_epoch.trial{j};
        end
        
        % Put data into EEGLab data structure
        clear EEG
        nChans = channels_number_cephalic;
        EEG.data = dat_temp; EEG.nbchan = nChans;
        EEG.srate = 512;
        EEG.chanlocs = channel_locations;
        EEG.trials = length(dataX_epoch.trial);
        EEG.pnts = size(dat_temp,2);
        EEG.xmin = -prestim;
        EEG.xmax = poststim;
        EEG.icaact = [];
        EEG.epoch = dataX_epoch.sampleinfo;
        EEG.setname = 'dataX_epoch';
        EEG.filename = ''; EEG.filepath = '';
        EEG.subject = ''; EEG.group = '';
        EEG.condition = ''; EEG.session = [];
        EEG.comments = 'preprocessed with fieldtrip';
        EEG.times = dataX_epoch.time{1};
        EEG.ref = []; EEG.event = [];
        EEG.icawinv = []; EEG.icasphere = [];
        EEG.icaweights = []; EEG.icaact = [];
        EEG.saved = 'no'; EEG.etc = [];
        EEG.specdata = []; EEG.icachansind = []; EEG.specicaact = [];
        
        % Select bad channel removal thresholds
        [~,idx1] = pop_rejchan(EEG,'elec',1:nChans,'threshold',3,...
            'norm','on','measure','kurt');
        [~,idx2] = pop_rejchan(EEG,'elec',1:nChans,'threshold',3,...
            'norm','on','measure','prob');
        [~,idx3] = pop_rejchan(EEG,'elec',1:nChans,'threshold',3,...
            'norm','on','measure','spec');
        badChans = unique([idx1,idx2,idx3]);
        
        % Spline interpolate bad channels
        if ~isempty(badChans)
            EEG = pop_interp(EEG,badChans,'spherical');
        end
        
        % Organise for Ft Struct
        dat_temp2 = cell(1,length(dataX_epoch.trial));
        for j = 1:length(dataX_epoch.trial)
            dat_temp2{j}(:,:) = double(EEG.data(:,:,j));
        end
        dataX_rm = dataX_epoch;
        dataX_rm.trial = dat_temp2;
        clear dat_temp2 dat_temp dataX_epoch
        
        %-------------------------------------
        %% Remove trials above STD THRESHOLD
        %-------------------------------------
        clear dat_temp
        for j = 1:length(dataX_rm.trial)
            dat_temp(:,:,j) = dataX_rm.trial{j};
        end
        
        %find maximum values in each epoch
        NumSTD = 4; %if max amplitude value is greater than 3 standard deviations of median(max value of every other trial),then remove it
        nChans = size(dataX_rm.label,1);
        thrMin = -prestim; %lower bound of timewindow to which STD-based artifcat rejection criteria will be applied
        thrMax = poststim; %changed by KW 1/31/19. This is the upperbound of the timewindow (0=stim onset) to which STD based artifact will be applied
        thrMin2 = -prestim;
        thrMax2 = poststim;
        
        % Define max value
        maxVals = [];
        maxVals = [maxVals;squeeze(max(max(abs(squeeze(dat_temp(:,:,:))))))];
        
        % Set threshold based on median absolute deviaiton from max values
        thr2 = round(median(maxVals) + NumSTD*median(abs(maxVals-median(maxVals))));
        
        %remove trials if they contain a value that exceeds thr2
        [trl_in, trl_out, cl_dat, nElec] = eegthresh(dat_temp,size(dat_temp,2),1:nChans,-thr2,thr2,[thrMin thrMax],thrMin2, thrMax2); %2/15/19 KW - moved this line so that trial rejection comes after baseline correction
        
        % Organise for Ft Struct
        dat_temp2 = cell(1,size(cl_dat,3));
        for j = 1:size(cl_dat,3)
            dat_temp2{j}(:,:) = cl_dat(:,:,j);
        end
        dataX_rm.trial = dat_temp2;
        clear dat_temp2 dat_temp dataX_epoch
        
        % Fix data structure to match new trials
        dataX_rm.trialinfo(trl_out,:) = [];
        dataX_rm.time(trl_out) = [];
        dataX_rm.sampleinfo(trl_out,:) = [];
        
        %% Downsample data
        cfg = [];
        cfg.resamplefs = 256;
        dataX_rm_ds = ft_resampledata(cfg,dataX_rm);
        clear dataX_rm
        
        %% Rereference data
        cfg = [];
        cfg.reref = 'yes';
        cfg.refchannel = 'all';
        cfg.refmethod = 'avg';
        dataX_clean_ReRef = ft_preprocessing(cfg,dataX_rm_ds);
        
        % Define neighbourhood structure
        elec = ft_read_sens([study_path,'/','Resources_Misc','/','BioSemi','_',num2str(channels_number_cephalic),'_','AB','.sfp'],'filetype','sfp');
        cfg=[];
        cfg.layout = ['biosemi',num2str(channels_number_cephalic),'.lay'];
        cfg.method = 'triangulation';
        cfg.elec = elec;
        neighbours = ft_prepare_neighbours(cfg);
        % ft_neighbourplot(cfg)
        
        %% Apply Spline CSD
        cfg = [];
        cfg.channel = 'EEG';
        cfg.neighbours = neighbours;
        cfg.layout = ['biosemi',num2str(channels_number_cephalic),'.lay'];
        cfg.trials = 'all';
        cfg.method = 'spline'; %'hjorth'; %'spline';
        cfg.lambda = 1e-05;
        dataX_CSD_spline = ft_scalpcurrentdensity(cfg,dataX_rm_ds);
        
        %% Append and index
        data_ReRef{block_idx} = dataX_clean_ReRef;
        data_CSD_spline{block_idx} = dataX_CSD_spline;
        data_rm_ds{block_idx} = dataX_rm_ds;
        %data_CSD_kayser{block_idx} = dataX_CSD_kayser;
        
    end
    
    %% Concatinate all recording blocks
    eeg_ReRef = ft_appenddata([],data_ReRef{:});
    eeg_CSD_spline = ft_appenddata([],data_CSD_spline{:});
    eeg_dataX_rm_ds = ft_appenddata([],data_rm_ds{:});
    %eeg_CSD_kayser = ft_appenddata([],data_CSD_kayser{:});
    clear data_ReRef data_CSD_spline data_CSD_hjorth
    
    %% Save ReRef
    % Verify Directory Exists and if Not Create It
    if exist([study_path,'/',study_name,'/','Recordings_preCl_ReRef_new','/',subjects{subject_idx},'/'],'dir') == 0
        mkdir([study_path,'/',study_name,'/','Recordings_preCl_ReRef_new','/',subjects{subject_idx},'/']);
    end
    % Save Figures and Data
    filename = [study_path,'/',study_name,'/','Recordings_preCl_ReRef_new','/',subjects{subject_idx},'/',...
        subjects{subject_idx},'_','ft_data_struct'];    filetype = '.mat';
    save([filename,filetype],'eeg_ReRef','-v7.3'); clear filename filetype
    clear eeg_ReRef
    
    %% Save CSD spline
    % Verify Directory Exists and if Not Create It
    if exist([study_path,'/',study_name,'/','Recordings_preCl_CSD_spline_new','/',subjects{subject_idx},'/'],'dir') == 0
        mkdir([study_path,'/',study_name,'/','Recordings_preCl_CSD_spline_new','/',subjects{subject_idx},'/']);
    end
    % Save Figures and Data
    filename = [study_path,'/',study_name,'/','Recordings_preCl_CSD_spline_new','/',subjects{subject_idx},'/',...
        subjects{subject_idx},'_','ft_data_struct'];    filetype = '.mat';
    save([filename,filetype],'eeg_CSD_spline','-v7.3'); clear filename filetype
    clear eeg_CSD_spline
    
    %% Save DS
    % Verify Directory Exists and if Not Create It
    if exist([study_path,'/',study_name,'/','Recordings_preCl_dataX_rm_ds_new','/',subjects{subject_idx},'/'],'dir') == 0
        mkdir([study_path,'/',study_name,'/','Recordings_preCl_dataX_rm_ds_new','/',subjects{subject_idx},'/']);
    end
    % Save Figures and Data
    filename = [study_path,'/',study_name,'/','Recordings_preCl_dataX_rm_ds_new','/',subjects{subject_idx},'/',...
        subjects{subject_idx},'_','ft_data_struct'];    filetype = '.mat';
    save([filename,filetype],'eeg_dataX_rm_ds','-v7.3'); clear filename filetype
    clear eeg_dataX_rm_ds
    
    %     %% Save CSD kayser
    %     % Verify Directory Exists and if Not Create It
    %     if exist([study_path,'/',study_name,'/','Recordings_preCl_CSD_kayser','/',subjects{subject_idx},'/'],'dir') == 0
    %         mkdir([study_path,'/',study_name,'/','Recordings_preCl_CSD_kayser','/',subjects{subject_idx},'/']);
    %     end
    %     % Save Figures and Data
    %     filename = [study_path,'/',study_name,'/','Recordings_preCl_CSD_kayser','/',subjects{subject_idx},'/',...
    %         subjects{subject_idx},'_','ft_data_struct'];    filetype = '.mat';
    %     save([filename,filetype],'eeg_CSD_kayser','-v7.3'); clear filename filetype
    %     clear eeg_CSD_kayser
    
    
end


