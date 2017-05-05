%% Set up root folders
mainFolder = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos';

% add utils to path
currentFolder = 'utils';
addpath(fullfile(mainFolder,currentFolder));

% add fastICA
addpath('D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\toolboxes\FastICA_2.5');

% move to tests folder
currentFolder = 'tests';
cd(fullfile(mainFolder,currentFolder));




%% load CartoXP
load(fullfile(mainFolder, 'data\cartoXP\handles.mat'))

information  = handles.DB;
data         = information.EG;
variables.fs = 1000;

%% POINT 116

% Aitor gets index 96 (not point) in page 41
idx     = find(information.npunto == 116);
ecgData = data(idx);

dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;

for idx=1:numel(dataIn.signalNames)
    fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
end

dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;
dataOut            = preprocessForICA (dataIn);


%% (v) Merge EGM and ECG: Use same leads as Aitor

lisftOfSignals = {'R1','R2','M1','M2','AVF','I','V1','V4','V6'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);


% survey data
for idx=1:numSignalsICA
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(dataForICA(idx, :), signalNamesICA{idx}, dataIn.fs);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    drawnow;
end


%% Run fastICA, sort by eKurtosis and plot signals
[signalICA, A, W] = fastica(dataForICA);

% get eKurtosis
kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);

icaSorted = signalICA(idxKurt, :);


for idx=1:numSignalsICA
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end

%%
clear all, close all, home