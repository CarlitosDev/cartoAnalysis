%% my code

addpath('D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\carlosUtils');
addpath('D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\toolboxes\FastICA_2.5');

% that loads 'handles'
load('D:\Work\MATLABPromotionsModel\Matlab PhD\TesisCarlos\3-CodigoIcaFA - Aitor\3-CodigoIcaFA\handles.mat')
fieldnames(handles)


information = handles.DB;
data        = information.EG;


%% Let's process ONLY ECG signals and do fastICA

% Aitor gets point 95 in page 40
idx     = find(information.npunto == 95);
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
lisftOfSignals = {'AVF','AVL','AVR','I','II','III','V1','V2','V3','V4','V5','V6'};
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