%% my code

addpath('C:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\carlosUtils');
addpath('C:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\toolboxes\FastICA_2.5');

% that loads 'handles'
load('E:\Matlab PhD\TesisCarlos\3-CodigoIcaFA - Aitor\3-CodigoIcaFA\handles.mat')
fieldnames(handles)


information = handles.DB;
data        = information.EG;
variables.fs = 1000;


%% Let's process ONLY ECG signals and do fastICA

% Aitor gets index 95 (not point) in page 40
idx     = find(information.npunto == 115);
ecgData = data(idx);

dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;
dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;
dataOut            = preprocessForICA (dataIn);


for idx=1:numel(dataIn.signalNames)
    fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
end

% Only ECG
lisftOfSignals = {'R1', 'R2', 'M1', 'M2', 'AVF', 'I', 'V1', 'V4', 'V6'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');

dataForICA = dataOut.signalData(idxSignals, :);

% plot original vs filtered
figure,
signalIdx   = 1;
plot(dataIn.signalData(signalIdx, :));
hold on
plot(dataOut.signalData(signalIdx, :), 'g');


getFloorPow2   = @(x) 2^(nextpow2(length(x))-1);
getQuarterPow2 = @(x) 2^(nextpow2(length(x))-2);
hammingWindowSize = getQuarterPow2(dataOut.signalData(signalIdx, :));


%% carlos add here
% [pxx, wxx] = getPowerSpectralDensity(dataOut.signalData(signalIdx, :), dataIn.fs, ...
%  'hammingWindow', hammingWindowSize, 'lengthFFT', 8192,  'numberOverlap', [], 'psdCutOff', 0.85);

[pxx,w] = pwelch(dataOut.signalData(signalIdx, :));

omega = w*dataIn.fs/pi;
figure, plot(omega, pxx)

% get 90% of the spectrum energy 
normPxx = pxx / sum(pxx);
cumNormPxx = cumsum(normPxx);
idx        = cumNormPxx < 0.95;

pxxCapped   = pxx(idx);
omegaCapped = omega(idx);
figure, plot(omegaCapped, pxxCapped)

%%

% run fastICA
[signalICA, A, W] = fastica (dataForICA);

% get eKurtosis
kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);

icaSorted = signalICA(idxKurt, :);




figure,
for idx = 1:numSignalsICA/2  
  subplot(numSignalsICA/2, 1, idx)
  plot(icaSorted(idx, :));
  title(sprintf('ICA k = %3.2f', eKurtValues(idx))) 
end

figure,
for idx = 1:floor(numSignalsICA/2)
  subplot(numSignalsICA/2, 1, idx)
  plot(icaSorted(idx+ceil(numSignalsICA/2), :));
  title(sprintf('ICA k = %3.2f', eKurtValues(idx+ceil(numSignalsICA/2))))
end


idx = 1;
[pxx,w] = pwelch(icaSorted(idx, :));

omega = w*dataIn.fs/pi;
% get 90% of the spectrum energy 
normPxx = pxx / sum(pxx);
cumNormPxx = cumsum(normPxx);
idx        = cumNormPxx < 0.95;

pxxCapped   = pxx(idx);
omegaCapped = omega(idx);
figure, plot(omegaCapped, pxxCapped)

idx = 1;
hist(icaSorted(idx, :));