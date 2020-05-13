function A = fltanlyteeg(eeg,Fs,cmpuse)
% Bandpass filter the signal using filtfilthd.
% For N EEG channels, where each channel is a different column in the eeg
% matrix, grouped by component
% The types of components included can also be explicitly specified in the cmpuse array:
%   1 = delta (0.2-4Hz), 2 = theta (4-8Hz), 3 = alpha (8-15Hz), 4 = beta
%   (15-30Hz)
% Nate Zuk (2017)

if nargin<3, cmpuse = 1:8; end % if not specified, use all components

% Preallocate arrays for storing the filtered EEG signals
if sum(cmpuse==1)>0, deltamag = NaN(size(eeg)); else; deltamag = []; end
if sum(cmpuse==2)>0, thetamag = NaN(size(eeg)); else; thetamag = []; end
if sum(cmpuse==3)>0, alphamag = NaN(size(eeg)); else; alphamag = []; end
if sum(cmpuse==4)>0, betamag = NaN(size(eeg)); else; betamag = []; end
if sum(cmpuse==5)>0, deltafns = NaN(size(eeg)); else; deltafns = []; end
if sum(cmpuse==6)>0, thetafns = NaN(size(eeg)); else; thetafns = []; end
if sum(cmpuse==7)>0, alphafns = NaN(size(eeg)); else; alphafns = []; end
if sum(cmpuse==8)>0, betafns = NaN(size(eeg)); else; betafns = []; end

for c = 1:size(eeg,2)
    % Create the bandpass filtered signals
    if sum(cmpuse==1)>0 | sum(cmpuse==5)>0 % delta
       flt = btrfiltereeg(eeg(:,c),Fs,'F3dB1',0.2,'F3dB2',4);
       hflt = hilbert(flt);
       if sum(cmpuse==1) % magnitude
           deltamag(:,c) = abs(hflt)';
       end
       if sum(cmpuse==5) % phase
           deltafns(:,c) = cos(angle(hflt));
       end
    end
    if sum(cmpuse==2)>0 | sum(cmpuse==6)>0 % theta
       flt = btrfiltereeg(eeg(:,c),Fs,'F3dB1',4,'F3dB2',8);
       hflt = hilbert(flt);
       if sum(cmpuse==2), % magnitude
           thetamag(:,c) = abs(hflt)';
       end
       if sum(cmpuse==6), % phase
           thetafns(:,c) = cos(angle(hflt));
       end
    end
    if sum(cmpuse==3)>0 | sum(cmpuse==7)>0 % alpha
       flt = btrfiltereeg(eeg(:,c),Fs,'F3dB1',8,'F3dB2',15);
       hflt = hilbert(flt);
       if sum(cmpuse==3) % magnitude
           alphamag(:,c) = abs(hflt)';
       end
       if sum(cmpuse==7) % phase
           alphafns(:,c) = cos(angle(hflt));
       end
    end
    if sum(cmpuse==4)>0 | sum(cmpuse==8)>0 % beta
       flt = btrfiltereeg(eeg(:,c),Fs,'F3dB1',15,'F3dB2',30);
       hflt = hilbert(flt);
       if sum(cmpuse==4) % magnitude
           betamag(:,c) = abs(hflt)';
       end
       if sum(cmpuse==8) % phase
           betafns(:,c) = cos(angle(hflt));
       end
    end
end

A = [deltamag thetamag alphamag betamag deltafns thetafns alphafns betafns];

unwrap