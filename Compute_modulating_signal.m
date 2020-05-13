% Summary:
% Script to extract envelope using Gammachirp filterbank

% Status:
% Under development

% Notes:
% cf. Crosse et al. 2017

% Author(s):
% Kevin Prinsloo

% Editor(s):
%

%% Prepare Workspace
%clearvars
%close all
%clc

%%  Script Config
% Dependencies
study_path = '/scratch/kprinslo/Chimera_Study';
path = '/scratch/kprinslo/Chimera_Study';

addpath(genpath('/scratch/kprinslo/Chimera_Study/Toolboxes/Envelope_extraction_tools'));
addpath '/scratch/kprinslo/Chimera_Study'
addpath '/scratch/kprinslo/Chimera_Study/Toolboxes/Envelope_extraction_tools';
addpath '/scratch/kprinslo/Chimera_Study/Scripts/';

% Define Envelope
envelope_type = 'GammaTone_256bands_Leagues_fix_fs'; % GammaToneFilterBank | HilbertEnvelope | GammaTone_256bands_TFS_Journey
trfEnvelope = 'none'; % 'trfChimaera'

% Verify Directory Exists and if Not Create It
if exist([path,'Stimuli','/',envelope_type,'/'],'dir') == 0
    mkdir([path,'Stimuli','/',envelope_type,'/']);
end

%% Envelope extraction
downsample_mod_sig = 'None'; % Downsampled | None
envelopes=cell(20,1);

%% Preallocation
audio_wav_gam1 = zeros(256,2646000);
audio_wav_gam1_TFS = zeros(256,2646000);

Fs = 44100;
Fmin = 80;	% lower frequency of filterbank in Hz
Fmax = .4*Fs;	% upper frequency of filterbank (.8 * Nyquist)
refilter = 1;
Nbands = 1;

% determine band cutoffs equally spaced on basilar membrane
Fco = equal_xbm_bands(Fmin, Fmax, Nbands);

audio_wav_holder = zeros(desired_sample_lengths,1);
Gammtone_filter_applicatione = 'Gram';

condition_idx_cluster = str2double(getenv('SLURM_ARRAY_TASK_ID'));  % this give's you back the job parameter from slurm (#SBATCH --array=1-16)
disp(condition_idx_cluster)


conditions = {'nb_0','nb_1','nb_4','nb_16'};
for condition_idx = 1:condition_idx_cluster
    condition_name = conditions{condition_idx};
    for trial_idx = 1:15, fprintf('Stimulus %d...\n',trial_idx),
        tic
        
        % Read in audio
        clear audio_wav       
        [wav_data_holder, ~] = audioread(['/scratch/kprinslo/Chimera_Study/Stimuli/Original_raw','/','Env_Leauges_fts_Journey_run_',num2str(trial_idx),'_',condition_name,'.wav']);
        
        % Select Journey trials from channel #2
        %leagues=wav_data_holder(:,1); 
        %leagues = resample(leagues,44100,fs);
        wav_data_holder = wav_data_holder(:,1);
        fs = 44100;
                   
        %% Compute multi-band chimeras - e1_fs2, e2_fs1 (env1_fts2_out)
        if strcmp(trfEnvelope,'trfChimaera')            
            [env_Leauges_fts_Journey_out, envjourney_ftsleauges_out] = multi_band_chimera(leagues, journey, Fco, Fs, refilter);
            tmp = hilbert(env_Leauges_fts_Journey_out);
            tfs_bands = cos(angle(tmp));
        elseif strcmp(trfEnvelope,'none')   
            tfs_bands = wav_data_holder;
        end
                
        %% GammacHirp filterbank
        GCparam.fs = fs;
        GCparam.NumCh = 256;
        GCparam.FRange = [80,3e3];
        GCparam.OutMidCrct = 'ELC';
        % GCparam.OutMidCrct = 'No';
        % GCparam.Ctrl = 'dyn';
        
        %% Filter below Nyquist frequency
        % Lowpass filter
        Fpass = 3.6e3;
        Fstop = 4.8e3;
        Fs=fs;
        Apass = 1;
        Astop = 60;
        h = fdesign.lowpass(Fpass,Fstop,Apass,Astop,fs); % h = fdesign.lowpass(Fpass,Fstop,Apass,Astop,fs);
        lpf1 = design(h,'cheby2','MatchExactly','stopband');
        clear Fpass Fstop Fs h;
        tfs_bands = filtfilthd(lpf1,tfs_bands);
        
        %% Apply Gammatone Filter band - Cochlear Processing
        tfs_bands=tfs_bands';
        audio_wav_gam = GCFBv210(tfs_bands,GCparam);
        clear audio_wav tfs_bands
        
        %% Calculate narrowband and broadband envelopes
        clear audio_wav_gam1
        for chn=1:size(audio_wav_gam,1)
            fprintf('Channel %d...\n',chn),
            audio_wav_gam1(chn,:)=abs(hilbert(audio_wav_gam(chn,:)));
        end        
        clear envelope1 envelope
        envelope1 = mean(audio_wav_gam1,1); clear audio_wav_gam1 audio_wav_gam

        %% Downsample                
        % Plot Original Data
        figure_1 = figure; plot(envelope1,'Color',[0 0.4470 0.7410]); axis tight; ylim([min(envelope1) max(envelope1)]);
        axis_1 = gca; axis_1.XColor = [0 0.4470 0.7410]; axis_1.YColor = [0 0.4470 0.7410];
        xlabel('Time (Sampling Periods)'); ylabel('Voltage (uV)');
        
        % Down sample signal
        envelope = resample(envelope1,128,fs);
        envelope = abs(envelope');    
        
        % Plot Downsampled Data
        axis_1_pos = axis_1.Position; axis_2 = axes('Position',axis_1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
        axis_2.XColor = 'r'; axis_2.YColor = 'r'; axis_2.ActivePositionProperty = 'outerposition';
        hold on; plot(envelope,'Parent',axis_2,'Color','r'); axis tight; ylim([min(envelope) max(envelope)]);
                     
        %envelope_norm = normaliseVector(envelope);        
        figure_2 = figure;
        plot(envelope);
               
        % Save Data
        filename = ([path,'Stimuli','/',envelope_type,'/','Env_Leauges_fts_Journey_run_',...
            num2str(trial_idx),'_',condition_name]);
        filetype = '.fig';
        savefig(figure_1,[filename,'_preFilter',filetype]); close(figure_1); clear figure_1
        savefig(figure_2,[filename,'_postFilter',filetype]); close(figure_2); clear figure_2 filetype
        filetype = '.mat';
        save([filename,filetype],'envelope','-v7.3');
        clear envelope filename filetype audio_wav fs spectrogram difference_samples desired_sample_lengths lpf1 h audio_wav_gam envelope 
        close all
        toc
    end
end


