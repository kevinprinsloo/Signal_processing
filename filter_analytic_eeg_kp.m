% Summary:
% Bandpass filter the signal using filtfilthd.

% Status:
% Under development

% Notes:
%  > For N EEG channels, where each channel is a different column in the eeg
%  > matrix, grouped by component
%  > The types of components included can also be explicitly specified in the temporal_fine_structure_index array:
%  > 1 = delta (0.2-4Hz) | 2 = theta (4-8Hz) | 3 = alpha (8-15Hz) | 4 = beta 15-30Hz)

% Author(s):
% Kevin Prinsloo

% Editor(s):
% 

function [eeg] = filter_analytic_eeg_kp(eeg_temp,Fs,filter_eeg_type_idx)
%function [analytic_output, hilbert_output] = filter_analytic_eeg(eeg,Fs,temporal_fine_structure_index)

for Channel_idx = 1:size(eeg_temp,2),fprintf('Channel %d',Channel_idx)
    fprintf('......')
    
    % Remove EEGlab - error with filtering
    rmpath(genpath('D:\Multisensory_Integration_Project\Toolboxes\eeglab_current\eeglab14_1_2b'))
    
    %% 
    if filter_eeg_type_idx == 1       
    BpFilt = designfilt('bandpassfir', ...
        'StopbandFrequency1',0.0001, ...
        'PassbandFrequency1',0.25, ...
        'PassbandFrequency2',1, ...
        'StopbandFrequency2',2, ...
        'DesignMethod','kaiserwin',...
        'PassbandRipple',0.01, ...
        'StopbandAttenuation1',50, ...
        'StopbandAttenuation2',50, ...
        'SampleRate',Fs);
    elseif filter_eeg_type_idx == 1  
    
    % Add EEGlab - error with filtering
    addpath(genpath('D:\Multisensory_Integration_Project\Toolboxes\eeglab_current\eeglab14_1_2b'));
    
    %fv1=fvtool(BpFilt)
    flteeg = filtfilt(BpFilt,eeg); %filtfilthd
    
    % Transform into hilbert complex
    eeg = hilbert(flteeg);
    
end

fprintf('Elapsed time: %0.1f minutes\n',toc/60)
end



