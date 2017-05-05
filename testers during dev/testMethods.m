%% Generate a 10Hz and 25Hz signal with noise and baseline drift
n   = 1:2500;
Fs  = 1e3;

% (A)
f1  = 10;
omega = 2.5*pi*f1/Fs;
x1 = 1.25*sin(omega * n);

% (B)
f2 = 28;
omega = 2*pi*f2/Fs;
x2 = 2.5*sin(omega * n);

% (C)
noise = 0.75*rand(1, length(n));

% (D)
omega = 2*pi*1.2/Fs;
baselineWander = 1.0*sin(omega * n);

% power line interference 50Hz
omega = 2*pi*50/Fs;
powerLineInterference = 0.25*sin(omega * n);

synSignal = x1+x2+noise+baselineWander+powerLineInterference;

figureRight,
plot(n, synSignal)
title('Input signal');
originalAxes = gca();

inputSignal  = synSignal;
samplingFreq = Fs;
lengthSignal = numel(inputSignal);

%% get the PSD for the current input signal
[inputPSD, inputPSDOmega] = ...
getPowerSpectralDensity(inputSignal, samplingFreq, ...
'windowType', 'rectangular'        , ...
'windowSize', length(inputSignal), ...
'lengthFFT' , 8192, 'numberOverlap', [], 'psdCutOff', 100);

figure,
plot(inputPSDOmega, inputPSD);
title('PSD Input signal');

%% Preprocess

dataIn.signalData = inputSignal;
dataIn.fs = samplingFreq;


% filter
f1HighPassCutOff = 1.0;
f2StopBandCutOff = 50 ;
f3LowPassCutOff  = 250;

dataForICA = preprocessForICAv2(dataIn, ...
    f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);
  
% detrend
%[dataForICA, baselineSignal] = removeBaselineWandering(dataForICA);
samplingPercentage = 3;
[dataForICA, baselineSignal] = removeBaselineWanderingV2(dataForICA, samplingPercentage);
  
hold(originalAxes, 'on');
plot(originalAxes, n, dataForICA, 'r')
title('Filtered Input signal');


figure,
plot(n, baselineWander)
hold on
plot(n, baselineSignal, 'r')
title('Baseline vs trend')


%%

[inputPSD, inputPSDOmega] = ...
getPowerSpectralDensity(dataForICA, samplingFreq, ...
'windowType', 'rectangular'        , ...
'windowSize', length(dataForICA), ...
'lengthFFT' , 8192, 'numberOverlap', [], 'psdCutOff', 100);

figure,
plot(inputPSDOmega, inputPSD);
title('PSD Filtered signal');

%%

[inputHatPSD, inputHatPSDOmega] = ...
    getPowerSpectralDensityBotteron(inputSignal, samplingFreq, ...
    'windowType'   , 'rectangular', ...
    'windowSize'   , lengthSignal/4, ...
    'lengthFFT'    , 8192, ...
    'numberOverlap', [], ...
    'psdCutOff'    , 60);

figure,
plot(inputHatPSDOmega, inputHatPSD);
[inputHatFreqPeak, inputHatFreqPeakIdx] = max(inputHatPSD);
maxFrequency = inputHatPSDOmega(inputHatFreqPeakIdx);

fprintf('C-Frequency Peak     %3.2f\n', maxFrequency);
%% Get funfamental freq
[fDom,indff,bwfromf0,z,periodogramInput,fz] = df_Ng(inputSignal, samplingFreq);
% deltaFreq = fz(2)-fz(1);
% ancho = 2*deltaFreq;
ancho = 1.8;%1Hz
[oi,ejexplot,ejeyplot] = organizationIndex(periodogramInput,fz,fDom,ancho);
figure, plot(ejexplot, ejeyplot)
fprintf('Frequency Peak     %3.2f\n', fDom);
fprintf('Organisation Index %3.2f\n', oi);