% Summary:
% Script to plto analyses from autotune MAMTA

% Status:
% Under Development

% Notes:
% n/a

% Author(s):
% Kevin Prinsloo

% Editor(s):
% 

%----------------------
%% Prepare Workspace
%----------------------
clearvars
close all
clc

%% Initialise Workspace
cd 'Y:\Kevin\Attention_Study\Scripts';
study_path = 'Y:\Kevin\Attention_Study\';
save_path = 'Y:\Kevin\Attention_Study\Stimuli\Save_stimuli\';

stim_path = 'Y:\Kevin\Attention_Study\111111111111111\SamHarris_c1\';

addpath 'Y:\Kevin\Attention_Study\Scripts\';
addpath 'Y:\Kevin\Attention_Study\ChimeraSoftware\';
addpath 'Y:\Kevin\Attention_Study\Toolboxes\HRTF-Individualization-master\HRTF-Individualization-master';
addpath 'Y:\Kevin\Toolboxes\CIPIC_hrtf_database\CIPIC_hrtf_database';
addpath 'Y:\Kevin\Multisensory_Integration_Project\Toolboxes\NoiseTools\NoiseTools';

%% Define workspace variables
frequency_sample = 44100; %48000;
trial_duration = 60;
samples_per_trial = round(frequency_sample*trial_duration); 
trial_numbers = 1;
Fs = 44100; %2881600;

clear aud_1 aud_2 tmp tmp2
% Define Obama trials
wav_files_obama = cell(trial_numbers);%zeros(trial_numbers,2881600,2);
aud_1 = cell(trial_numbers,1);

for itrials = 1   %:trial_numbers
    
    % Load wav files - data [2881600 x 1]
    [wav_data_holder, fs] = audioread([stim_path,'SamHarris_c1_',num2str(itrials),'.mp4']);
    
end


% Define Bush trials
clear tmp tmp2
wav_files_bush = zeros(size(aud_1{itrials}));
aud_2 = cell(itrials,1);
for itrials = 3%:trial_numbers
    [wav_data_bush, FS_k] = audioread([study_path,'Stimuli','\','Audio_Files_Bush','\','bush_a_',num2str(itrials),'.wav']);
    tmp = length(aud_1{itrials}(:,:));
    for ichan = 1:2
        tmp2 = (length(wav_data_bush(:,ichan)))-tmp;
        wav_files_bush(:,ichan) = wav_data_bush(1:end-tmp2,ichan);
    end
    aud_2{itrials} = wav_files_bush; clear wav_data
end

% Do a quick check of the audio files loaded
gong = audioplayer(aud_2{itrials}, FS_k);
play(gong);
stop(gong);

% This is set if using the toolbox
%Nbands = 3;
do_play = 0;
do_plot = 1;

% Create Auditory Chimera with the chosen trials
Stero_channel_chimera = cell(trial_numbers);
nbands_list = 1;
for itrial = 3    
    Nbands = 1;
    
    Fmin = 80;	% lower frequency of filterbank in Hz
    Fmax = .4*Fs;	% upper frequency of filterbank (.8 * Nyquist)
    refilter = 1;
    
    % determine band cutoffs equally spaced on basilar membrane
    Fco = equal_xbm_bands(Fmin, Fmax, Nbands);
    
    % compute multi-band chimeras - e1_fs2, e2_fs1 (env1_fts2_out)
    clear audio_1 audio_2
    audio_1 = aud_1{itrial}; audio_2 = aud_2{itrial};    
    [env1_fts2_out, env2_fts1_out] = multi_band_chimera(aud_1{itrial}, aud_2{itrial}, Fco, Fs, refilter);
           
    % Define channels
    %>> By convention, the left channel is the first column and the right channel is the second column.
    %tmp = zeros(length(env2_fts1_out),1);
    %Stero_channel_left = [env2_fts1_out tmp];
    %Stero_channel_chimera_tmp = [env2_fts1_out wav_files_bush];
    
    % normalize and save
    %env1_fts2 = env1_fts2_out./max(abs(env1_fts2_out(:)));
    Stero_channel_chimera{itrial} = env2_fts1_out./max(abs(env2_fts1_out));
    
    %clear Stero_channel_chimera
    %for k = 1:2
    %    Stero_channel_chimera(:,k) = Stero_channel_chimera_tmp(:,k)./max(abs(Stero_channel_chimera_tmp(:,k)));
    %end
    
%     filetype = '.wav';
%     %filename = ([save_path,'\','env1_fts2','_',num2str(Nbands)]);
%     %audiowrite([filename,filetype],env1_fts2,frequency_sample);
%     filename = ([save_path,'\','Chimera_env2_fts1','_',num2str(Nbands)]);
%     audiowrite([filename,filetype], Stero_channel_chimera,frequency_sample);
    
end

%% Testing - play chimera
Stero_channel_chimera_test = Stero_channel_chimera{itrials};
gong = audioplayer(Stero_channel_chimera_test, FS_k);
play(gong);
stop(gong);

% Save Chimera wav file
filetype = '.wav';
filename = ([save_path,'/','Chimera_only','/','Chimera_env2_fts1','_','bnd','_',num2str(Nbands),'_','vid_' num2str(itrials)]);
audiowrite([filename,filetype], env2_fts1_out,frequency_sample);

% Define two audio channels
clear audioChan_1 audioChan_2
audioChan_1 = (Stero_channel_chimera{itrial}); 
audioChan_2 = (aud_2{itrial}); 

dir='Y:\Kevin\Attention_Study\';
Audiopath = [dir 'Audiobooks/Stimuli/Save_stimuli/'];    % plain audiofiles
%HRIR_path = [dir 'Toolboxes/CIPIC_hrtf_database/CIPIC_hrtf_database/Kemar/Small_pinna/'];  % standard KEMAR mannequin HRTF
HRIR_path = [dir 'Toolboxes\CIPIC_hrtf_database\CIPIC_hrtf_database\special_kemar_hrir\kemar_horizontal\'];
Savepath = [dir 'Stimuli/HRTF_applied/']; % small pinna folder for subject_165 or large pinna for subject 021

% Load in the HRIR 
% corresponds to standard KEMAR small pinna
% and the HRIR is 200 samples in length corresponding to ~4.5 ms (see CIPIC_HRTF_Database.pdf)
tmp = load([HRIR_path 'small_pinna_final.mat'],'left','right');  % corresponds to standard KEMAR large pinna
%>> left: [200×72 double] | right: [200×72 double]
hrir_left_30 = left(:,7); hrir_right_30 = right(:,7);
hrir_left_zero = left(:,1); hrir_right_zero = right(:,1);

%----------------
%% Chimera
%---------------
clear rms_aud Conv_left Conv_right HRTF_aud_holder HRTF_aud_holder
%divide by rms and then by max before writing
rms_aud = [(audioChan_1(:,1)./nt_rms(audioChan_1(:,1))),(audioChan_1(:,2)./nt_rms(audioChan_1(:,2)))]; %normalise and concatenate
%convolve the left HRIR with channel 1 and right HRIR with channel 2 of the audio
Conv_left = conv(rms_aud(:,1),hrir_left_zero,'same'); % Channel 1 - left - obama
Conv_right = conv(rms_aud(:,2),hrir_right_zero,'same'); % Channel 2 - right
HRTF_aud_holder = [Conv_left Conv_right];
%divide by max to prevent clipping when writing to file
HRTF_aud_chimera = HRTF_aud_holder./max(abs(HRTF_aud_holder));
clear rms_aud Conv_left Conv_right HRTF_aud_holder HRTF_aud_holder

%----------------
%% Attending
%---------------
clear rms_aud Conv_left Conv_right HRTF_aud_holder HRTF_aud_holder
%divide by rms and then by max before writing
rms_aud = [(audioChan_2(:,1)./nt_rms(audioChan_2(:,1))),(audioChan_2(:,2)./nt_rms(audioChan_2(:,2)))]; %normalise and concatenate
%convolve the left HRIR with channel 1 and right HRIR with channel 2 of the audio
Conv_left = conv(rms_aud(:,1),hrir_left_30,'same'); % Channel 1 - left - obama
Conv_right = conv(rms_aud(:,2),hrir_right_30,'same'); % Channel 2 - right
HRTF_aud_holder = [Conv_left Conv_right];
%divide by max to prevent clipping when writing to file
HRTF_aud_attending = HRTF_aud_holder./max(abs(HRTF_aud_holder));
clear rms_aud Conv_left Conv_right HRTF_aud_holder HRTF_aud_holder

%-------------
%% Combining
%-------------
HRTF_combined = zeros(size(HRTF_aud_chimera));
for ichan = 1:2    
    HRTF_combined(:,ichan) = HRTF_aud_chimera(:,ichan)+HRTF_aud_attending(:,ichan);
end
HRTF_combined = HRTF_combined./max(abs(HRTF_combined));

filetype = '.wav';
filename = ([save_path,'\','Chimera_env2_fts1','_','bnd','_',num2str(Nbands),'_','vid_' num2str(itrials)]);
audiowrite([filename,filetype], HRTF_combined,frequency_sample);

% Test play HRTF Chimera
gong = audioplayer(HRTF_combined, FS_k);
play(gong);
stop(gong);

orig = HRTF_combined;

% Save Chimera
%write the audio file **check saving name matches the HRTF read in and the folder name
audiowrite([Savepath audio_channel{iAudio} '_Subject_021_30_0_' num2str(itrial) '.wav'],HRTF_aud,Fs);



















%%>> NOTES
%>> In Adobe After Effects CC 2018 - when rendering the audio+video avi -
%>> ensure to chnage best settings to current settings
%>> AV output will be ~1.25 GB (if not it will save up to 10 GB per avi.)

%-------------------------------------------
%% Head-related transfer functions (HRTFs)
%-------------------------------------------

% Apply HRTF to present source at 30 degrees azimuth
% uses HRTF for mannequin from the CIPIC database

dir='Y:\Kevin\Attention_Study\';

Audiopath = [dir 'Audiobooks/Stimuli/Save_stimuli/'];    % plain audiofiles
HRIR_path = [dir 'Toolboxes/CIPIC_hrtf_database/CIPIC_hrtf_database/Kemar/Small_pinna/'];  % standard KEMAR mannequin HRTF
Savepath = [dir 'Stimuli/HRTF_applied/']; % small pinna folder for subject_165 or large pinna for subject 021

% From CIPIC Database
% HRIR_path = [dir 'CIPIC Database/special_kemar_hrir/kemar_horizontal/'];  % special KEMAR mannequin HRTF

audio_channel = {'Attend','Unattend'};               %
nTrials = 18;

% Load in the HRIR 
load([HRIR_path 'hrir_final.mat'],'hrir_r','hrir_l');  % corresponds to standard KEMAR large pinna
% load([HRIR_path 'Subject_165_30_0.mat'],'hrir_right','hrir_left');  % corresponds to standard KEMAR small pinna
% and the HRIR is 200 samples in length corresponding to ~4.5 ms (see CIPIC_HRTF_Database.pdf)

%% Apply HRTF
for itrial = 1%:nTrials
    %for iAudio = 1:numel(audio_channel)

        % Load in the audio file to be filtered - this is the audiofile which we will extract the envelope from
        aud_1 = ([Audiopath audio_channel{itrial} '_run_' num2str(itrial) '.wav']);
        [y_aud_1,Fs] = audioread(aud_1);
        
        aud_2 = ([Audiopath audio_channel{itrial} '_run_' num2str(itrial) '.wav']);
        [y_aud_2,Fs] = audioread(aud_2);
                
        %convolve the left HRIR with channel 1 and right HRIR with channel 2 of the audio
        new_left = conv(y(:,1),hrir_left,'same');
        new_right = conv(y(:,2),hrir_right,'same');
        
        %divide by rms and then by max before writing
        HRTF_aud = [(new_left./nt_rms(new_left)),(new_right./nt_rms(new_right))]; %normalise and concatenate 

        %divide by max to prevent clipping when writing to file
        HRTF_aud = HRTF_aud./max(abs(HRTF_aud)); 
        %>> Apply log normilisation here? - as frequency distribution is
        %>> non-linear
       
        %write the audio file **check saving name matches the HRTF read in and the folder name
        audiowrite([Savepath audio_channel{iAudio} '_Subject_021_30_0_' num2str(itrial) '.wav'],HRTF_aud,Fs);

    %end 
end

 powspctrm = abs(data.fourierspctrm).^2;
 
 logspctrm = 10*log10(powspctrm./repmat(nanamean(powspctrm,4),[1 1 1 numel(data.time)]));
 ctmp.logspectrm = sequeeze(mean(logspctrm(:,:,:,:)));


%---------------
%% Graveyrad
%---------------

videoFWriter = vision.VideoFileWriter(fullfile(shotPath, 'avclip', ...
    [num2str(shotNum) '_av_clip.avi']), ...
    'AudioInputPort', true, ...
    'FrameRate',  annot.video.frame_rate);

frames = 2881600;
for itrial = 1:length(frames)
    %fprintf('Frame: %d/%d\n', i, length(frames));
    step(videoFWriter, frames{itrial}, subSignal(numRep*(itrial-1)+1:numRep*itrial,:));
end
release(videoFWriter);

%% NOTES

% Vid1_tfs2_Env1 - latsest version: Don't want to unattend AV stimuli to
% distract from understanding the words 
% What the stim to be distracting but not comprehendable 







