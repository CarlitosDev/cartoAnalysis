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

%% POINT 84 (medium atrium)

% Let's try medium atria
idx     = find(information.npunto == 84);
ecgData = data(idx);

dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;

for idx=1:numel(dataIn.signalNames)
    fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
end

dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;
dataOut            = preprocessForICA (dataIn);


%% (i.1) Only monopolar catheter's 
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


%% (i.2) Only monopolar catheter's named M
lisftOfSignals = {'M1','M2','M3','M4'};
numSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);

vaioFlag = true;
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
for idx=1:numSignalsICA
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end



% stupid test - compare a lead with the extracted ica from EGM

lisftOfSignals = {'AVF'};
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataToCompare     = dataOut.signalData(idxSignals, :);
dataToCompareName = dataIn.signalNames(idxSignals);

normToSum = @(x) x/sum(x);
normToMax = @(x) x/abs(max(x));
normBoth  = @(x) normToSum(normToMax(x));

idx = 4;

figure,
plot(normBoth(icaSorted(idx, :)));
hold on
plot(normBoth(dataToCompare), 'g');

outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
legend(outSignalName, dataToCompareName{1});

idx = 1;

figure,
plot(normBoth(icaSorted(idx, :)));
hold on
plot(normBoth(dataToCompare), 'g');

outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
legend(outSignalName, dataToCompareName{1});


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

%% (iii.a) Only ECG

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
for idx=1:numel(eKurtValues);
    outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
    [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA(icaSorted(idx, :), outSignalName, dataIn.fs);
    
    if vaioFlag, set(gcf,'units','normalized','outerposition',[0 0 1 1]),end;
end



%% (iii.b) Only ECG where AA is meant to be maximal

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


%% (iv) Merge EGM and ECG: This is the quid

lisftOfSignals = {'M1','M2','M3','M4', 'AVF','AVL','III','V1'};
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