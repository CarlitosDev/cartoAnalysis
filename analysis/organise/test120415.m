mainFolder = 'E:\Matlab PhD\matlab code';

% add utils to path
currentFolder = 'utils';
addpath(fullfile(mainFolder,currentFolder));

% add fastICA
addpath('E:\Matlab PhD\matlab code\toolboxes\FastICA_2.5');

% move to tests folder
currentFolder = 'tests';
cd(fullfile(mainFolder,currentFolder));




%% load CartoXP
load('E:\Matlab PhD\matlab code\data\cartoXP\handles.mat')

information  = handles.DB;
data         = information.EG;
variables.fs = 1000;

%% POINT 48 (medium atrium)

% Aitor gets index 95 (not point) in page 40
% Let's try medium atria
idx     = find(information.npunto == 48);
ecgData = data(idx);

dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;

for idx=1:numel(dataIn.signalNames)
    fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
end

dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;
dataOut            = preprocessForICA (dataIn);


%% (i) Only monopolar catheter's 
lisftOfSignals = {'M1','M2','M3','M4','R1','R2'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);


% survey data
for idx=1:numSignalsICA
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(dataForICA(idx, :), signalNamesICA{idx}, dataIn.fs);
end


% run fastICA
[signalICA, A, W] = fastica(dataForICA);

% get eKurtosis
kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);

icaSorted = signalICA(idxKurt, :);

% have a look to the output
vaioFlag = 1;
for idx=1:numSignalsICA
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end


%% (ii) Bipolar catheter's

lisftOfSignals = {'M1_MINUS_M2','M3_MINUS_M4','R1_MINUS_R2'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);


% survey data
for idx=1:numSignalsICA
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(dataForICA(idx, :), signalNamesICA{idx}, dataIn.fs);
end


% run fastICA
[signalICA, A, W] = fastica(dataForICA);

% get eKurtosis
kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);

icaSorted = signalICA(idxKurt, :);

% have a look to the output
vaioFlag = 1;
for idx=1:numSignalsICA
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end

%% Only ECG

lisftOfSignals = {'AVF','AVL','AVR','I','II','III','V1','V2','V3','V4','V5','V6'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);


% survey data
for idx=1:numSignalsICA
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(dataForICA(idx, :), signalNamesICA{idx}, dataIn.fs);
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end


% run fastICA
[signalICA, A, W] = fastica(dataForICA);

% get eKurtosis
kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);

icaSorted = signalICA(idxKurt, :);

% have a look to the output
vaioFlag = 1;
for idx=1:numSignalsICA
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end



%% Only ECG where AA is meant to be maximal

lisftOfSignals = {'AVF','AVL','III','V1'};
%lisftOfSignals = {'AVF','AVL','II','III','V1'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);


% survey data
for idx=1:numSignalsICA
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(dataForICA(idx, :), signalNamesICA{idx}, dataIn.fs);
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end


% run fastICA
[signalICA, A, W] = fastica(dataForICA);

% get eKurtosis
kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);

icaSorted = signalICA(idxKurt, :);

% have a look to the output
vaioFlag = 1;
for idx=1:numSignalsICA
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end


%% Merge EGM and ECG

lisftOfSignals = {'M1','M2','M3','M4', 'AVF','AVL','III','V1'};
%lisftOfSignals = {'AVF','AVL','II','III','V1'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);


% survey data
for idx=1:numSignalsICA
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(dataForICA(idx, :), signalNamesICA{idx}, dataIn.fs);
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end


% run fastICA
[signalICA, A, W] = fastica(dataForICA);

% get eKurtosis
kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);

icaSorted = signalICA(idxKurt, :);

% have a look to the output
vaioFlag = 1;
for idx=1:numSignalsICA
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end