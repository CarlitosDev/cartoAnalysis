% generateMainFreqMapsBotteron
% Use Botteron-Smith to calculate the main frequency map for a given EGM


%% Load data
%path2CartoFile     = fullfile('data','Auriculas','#2', 'Carto', 'EG.mat');
%path2CartoFile     = fullfile('data','Auriculas','#4','Carto1','EG.mat');
% path2CartoFile     = fullfile('data','Auriculas','#4', 'Carto2', 'EG.mat');
path2CartoFile     = fullfile('data','Auriculas','#5', 'Carto', '1-Map', 'EG.mat');
fullPath2CartoFile = fullfile(rootFolder(), path2CartoFile); 


[data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile);


%% Set the input signal and the points to analyse

% So far just use one signal
 
lisftOfSignals = {'M1_MINUS_M2'};
variables.fs   = 1e3;
listOfPoints   = 'all';


%% Set the preprocess parameters
% Schilling (p145) sets filters to (30, 150)
% Also cancel out 50Hz
f1HighPassCutOff = 2.5;
f2StopBandCutOff = 50 ; % power line interference (notch)
f3LowPassCutOff  = 150;

%%

allPoints = unique(pointsInfo);  
numPoints = numel(allPoints);

useSignalFlag = char(1, numPoints);
scaleFrom0To1 = @(x) (x-min(x))/(max(x)-min(x));


%% Allocate output data

mainFrequency = zeros(1, numPoints);

%% loop through the adquisition points

%for idxPoint = 26
for idxPoint=1:numPoints
    
currentResults = [];
currentPoint   = allPoints(idxPoint);
  
idx      = find(pointsInfo == currentPoint);
ecgData  = data(idx);
    

% fprintf('\nProcessing point %d\n', currentPoint);
   

%% Check all signals are present and preprocess input data.
dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;
dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;

[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);
assert(all(signalPresent), 'Missing signals');

selDataIn.signalData  = dataIn.signalData(idxSignals, :);
selDataIn.signalNames = dataIn.signalNames(idxSignals)  ;
selDataIn.fs          = dataIn.fs;


preprocessData = preprocessForICAv2(selDataIn, ...
    f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);


% Overlay Botteron's spectrum
[powSpectDens, psdOmega] = ...
    getPowerSpectralDensityBotteron(preprocessData, selDataIn.fs, ...
    'windowType', 'rectangular'          , ...
    'windowSize', length(preprocessData) , ...
    'lengthFFT' , 8192, 'numberOverlap'  , [], ...
    'psdCutOff' , 40);

[~, freqPeakIdx] = max(powSpectDens);

mainFrequency(idxPoint)  = psdOmega(freqPeakIdx);

fprintf('Point %d has peak frequency of %3.2f Hz\n', currentPoint, psdOmega(freqPeakIdx));

end


plotMeshAndEGMvalues( meshData, cartoData, mainFrequency );

figTitle = strrep(path2CartoFile, '/', ' ');
title(figTitle);

maxFreq  = max(mainFrequency);
minFreq  = min(mainFrequency);
xlabel(sprintf('Blue %3.2f Hz Red %3.2f Hz', minFreq, maxFreq));

