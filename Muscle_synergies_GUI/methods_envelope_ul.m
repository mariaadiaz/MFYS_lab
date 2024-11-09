% function [envelope, outliers] = methods_envelope(emg_signal,fs,method,spikes);
%
% This function calculates the envelope of an emg signal
%
% Inputs:
%         emg_signal - One emg signal    
%         fs         - Sample frequency  
%         method     - Two methods available:
%                      1) Usual: The most common process to calculate 
%                         the envelope from an emg. You can find more 
%                         details in the following package developed 
%                         for Python: pyemgpipeline
%                         pyemgpipeline link: A Python package for EMG processing
%                      2) Literature: Other studies have used that method. 
%                         The reference of the studies are below.
%                         - Jackson, R. W., & Collins, S. H. (2019). 
%                           Heuristic-based ankle exoskeleton control for 
%                           co-adaptive assistance of human locomotion. 
%                           IEEE trans. neural sys. and rehab. engineering
%        spikes     -  Select 'nospikes' if you want to delete the spikes 
%                      from your emg signal. Otherwise input 'spikes' or '0'.
% Outputs:
%         envelope   - QRS complex merge in one signal
%         outliers   -
%
% Maria Alejandra Diaz
% VUB, 2022
% ma.diaz@vub.edu
%
function [envelope, outliers, freq, power] = methods_envelope_ul(emg_signal,fs,method,spikes)

fnyq = round(fs)/2;

if strcmp(method,'literature') == 1
    signal = emg_signal - mean(emg_signal);    
    % Filter from Scano 2018, Robotic assistance for UL
    [b,a] = butter(5,50/fnyq,'high');
    hp_signal = filtfilt(b,a,signal);
%     [b_high,a_high] = butter(5,0.008/fnyq,'high'); % Delete the DC offset 0HZ
    [b_low,a_low] = butter(5,10/fnyq,'low'); % 10 or 6
    rect_signal = abs(hp_signal);
%     envelope2 = filtfilt(b_high,a_high,rect_signal);
    envelope = filtfilt(b_low,a_low,rect_signal);

    % Smooth the signal using a moving average filter
    windowWidth = round(0.2*fs); % 200 ms approx
    kernel = ones(windowWidth,1) / windowWidth;
    envelope = filter(kernel, 1, envelope);    
    
elseif strcmp(method,'usual') == 1
    
    signal = emg_signal - mean(emg_signal);
    [b,a] = butter(4,[30 450]/fnyq,'bandpass');
    bp_signal = filtfilt(b,a,signal); 
    rect_signal = abs(bp_signal);
    [b,a] = butter(4,20/fnyq,'low');
    envelope = filtfilt(b,a,rect_signal);
    
    % Smooth the signal using a moving average filter
    windowWidth = round(0.2*fs); % 200 ms approx
    kernel = ones(windowWidth,1) / windowWidth;
    envelope = filter(kernel, 1, envelope);
end

%------ Analysis in the frequency domain
win_length= 1; % seconds; freq. resolution 2/Tw=0.4 Hz
nrsamp=floor(win_length*fs);
overlap = 0.10;
nfft=2^(2+nextpow2(length(signal)));
[Pxx,fxx] = pwelch(signal,hann(nrsamp),[],[],fs);
ave_freq = sum(fxx.*Pxx)/sum(Pxx);
[freq,power] = meanfreq(Pxx,fxx);

% Delete spikes from the envelope if necessary
if strcmp(spikes,'nospikes') == 1
    [envelope,outliers] = deletespikes(envelope);
else
    outliers = [];
end

end