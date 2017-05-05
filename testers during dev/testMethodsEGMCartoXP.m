

%idxPoint=217;
%idxPoint=20;
idxPoint=10;


listOfSelectedIcs = {1};

%listOfSelectedIcs = {[1,2]};
%listOfSelectedIcs = {1, 2, 3, [1,2], [1,2,3]};

testType       = 'BotteronSmith';
%srcDescription = fullfile('#2', 'Carto');
% srcDescription = fullfile('#4', 'Carto1');
% srcDescription = fullfile('#4', 'Carto2');
srcDescription = fullfile('#5', 'Carto', '1-Map');

path2CartoFile     = fullfile('data','Auriculas',srcDescription, 'EG.mat');
fullPath2CartoFile = fullfile(rootFolder(), path2CartoFile); 


numFlyBys = 1;
doRemoveBaselineWandering = false;

% Do whitening and filtering if requested    
doPreprocess = true;

% define input signals
lisftOfSignals = {'M1_MINUS_M2'};
listOfPoints   = 'all';
        
    
%%   Add paths for root, results and utils folders. Also set fastICA path
%   and the mat file where the data is stored.

commaString = @(listNumbers) strcat(sprintf('%d,', listNumbers(1:end-1)), sprintf('%d', listNumbers(end)));   
signalsStr  = @(listSignals) [sprintf('%s,', listSignals{1:end-1}), sprintf('%s', listSignals{end})];
floatStr    = @(listNumbers) strcat(sprintf('%3.3f,', listNumbers(1:end-1)), sprintf('%3.3f', listNumbers(end)));   


resultsFolder = fullfile(rootFolder(), 'results2Points', testType);
if ~isdir(resultsFolder)
    mkdir(resultsFolder);
end


fileName = fullfile(resultsFolder, ['log', datestr(now), '.csv']);
fid = 1;
%fid = fopen(fileName , 'a');
    
%% Set up configuration

%   Set the sampling frequency, the list of input signals for fastICA and
%   the points to extract the data from. The latter can be just a list of
%   points or the string 'all' meaning that the process will be repeated
%   through all the data.

    variables.fs = 1000;   
    samplingFreq = 1000;   
        

%% load CartoXP

[data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile);



%% where to apply fastICA.

if isa(listOfPoints, 'char') && strcmpi(listOfPoints, 'all')
    allPoints = unique(pointsInfo);  
    % let's assume numeric otherwise    
else
    % check the requested point belongs to the list. Otherwise, just wipe
    % it out.  
    currentPoints = unique(pointsInfo);
    presenPoints  = ismember(listOfPoints, currentPoints);
    allPoints     = listOfPoints(presenPoints);   
    if listOfPoints(~presenPoints)
        fprintf('Warning: Point %d not found\n', listOfPoints(~presenPoints));
    end
end

numPoints      = numel(allPoints);
numSelectedICs = numel(listOfSelectedIcs);

botteronFrequency = zeros(1, numPoints);
cartoPointAnalysis = cell(1, numPoints);

% Populate the header
fprintf(fid, '%s\t%s\t%s\t%s\ts%s\t%s\n', ...
    'Point', 'Input signals', 'input eKurtosis', ...
    'number of ICs', 'output eKurtosis', 'output peak Freq');

currentResults = [];
currentPoint   = allPoints(idxPoint);
  
idx      = find(pointsInfo == currentPoint);
ecgData  = data(idx);
    

fprintf('\nProcessing point %d\n', currentPoint);
    


% Check all signals are present and then preprocess input data.
dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;
dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;


[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);
assert(all(signalPresent), 'Missing signals');

selDataIn.signalData  = dataIn.signalData(idxSignals, :);
selDataIn.signalNames = dataIn.signalNames(idxSignals);
selDataIn.fs          = dataIn.fs;

% Preprocess
% Schilling (p145) sets filters to (30, 150)
[numSignals, lengthSignal] = size(selDataIn.signalData);

%%
inputSignal = selDataIn.signalData;

% (A)

lengthSignal = numel(inputSignal);
n = 1:lengthSignal;

figureRight,
plot(n, inputSignal)
originalAxes = gca();



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

% filter
f1HighPassCutOff = 1.0;
f2StopBandCutOff = 50 ;
f3LowPassCutOff  = 250;

dataForICA = preprocessForICAv2(selDataIn, ...
    f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);
  
% detrend
%[dataForICA, baselineSignal] = removeBaselineWandering(dataForICA);
samplingPercentage = 3;
[dataForICA, baselineSignal] = removeBaselineWanderingV2(dataForICA, samplingPercentage);
  
hold(originalAxes, 'on');
plot(originalAxes, n, dataForICA, 'r')
title('Filtered Input signal');



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