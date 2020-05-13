% Summary:
% Script to plot single-sided amplitude spectrum of data

% Status:
% Complete

% Notes:
% n/a

% Author(s):
% Kebin Prinsloo

function [f,pxx] = calcPSD(data,frequency_resolution,eeg_sampling_rate_downsampled_Hz)
window_length_samples = round(eeg_sampling_rate_downsampled_Hz/frequency_resolution);
if window_length_samples > length(data)
    window_length_samples = length(data);
end
noverlap = round(window_length_samples/2); % 50% overlap
nfft = 2^nextpow2(window_length_samples);
%[pxx,f] = pwelch(data,window_length_samples,noverlap,nfft,sampling_rate_Hz);
[pxx,f] = pwelch(data, hanning(512),0, nfft, eeg_sampling_rate_downsampled_Hz);
end


