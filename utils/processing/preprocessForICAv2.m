function dataOut = preprocessForICAv2(dataIn, ...
    f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff)

%% PREPROCESSFORICAV2 prepare EGM data for analysis or processing.
% 
%   After a nice talk with Doctor Charlie Butcher (Harefield-London), let's
%   be less severe with the filtering step.
%
% The preprocess consists of a high pass filter (1.0 Hz), followed by a
% notch filter to remove power line interference and low pass filter with a
% 250 Hz cut-off.
%
%   
%
% Carlos Aguilar - 15th November 2k15
%
% See also preprocessForICA.m


inputSignal  = dataIn.signalData;
samplingFreq = dataIn.fs;

% set preprocessing parameters.
% f1HighPassCutOff = 1.0;
% f2StopBandCutOff = 50 ;
% f3LowPassCutOff  = 250;

% Remove low frequency components. Let's set the order to 512(!!) to achieve
% at least -3dB

dataSignalProc = filterFIR(inputSignal, samplingFreq, ...
                'filterOrder'     , 512             , ...
                'filterType'      ,'high'           , ...
                'filterFrequency' , f1HighPassCutOff, ...
                'showFilterDesign', false);


% Use a Butterworth desing for cancelling electrical noise.
buttFilter = designfilt('bandstopiir'        , ...
                        'FilterOrder', 2     , ...
                        'HalfPowerFrequency1', f2StopBandCutOff-1,...
                        'HalfPowerFrequency2', f2StopBandCutOff+1,...
                        'DesignMethod'       , 'butter', ...
                        'SampleRate'         , samplingFreq);

% Gain zero-phase by using filtfilt.
dataSignalProc = filtfilt(buttFilter, dataSignalProc');

% 
dataSignalProc = filterFIR(dataSignalProc', samplingFreq, ...
                'filterOrder'     , 512                , ...
                'filterType'      ,'low'               , ...
                'filterFrequency' , f3LowPassCutOff    , ...
                'showFilterDesign', false);
            
% remove the mean
dataOut = ...
    bsxfun(@minus, dataSignalProc, mean(dataSignalProc,2));