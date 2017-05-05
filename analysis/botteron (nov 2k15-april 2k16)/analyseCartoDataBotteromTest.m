%% analyseCartoDataBotterom
% 
% Run Botteron-Smith analysis over a collection of carto points (EGM M1-M2)
% and bound the underlying geometry with the main frequency found in the
% analysis. Interact with the results by clicking on them to see the source
% data (both time and frequency).
% 
% Carlos Aguilar - February 2k16
%%

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




for idxPoint=1:numPoints

currentResults = [];
currentPoint   = allPoints(idxPoint);
  
idx      = find(pointsInfo == currentPoint);
ecgData  = data(idx);
    

fprintf('\nProcessing point %d\n', currentPoint);
    


%% Check all signals are present and then preprocess input data.
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

f1HighPassCutOff = 1.0;
f2StopBandCutOff = 50 ;
f3LowPassCutOff  = 250;
dataForICA = preprocessForICAv2(selDataIn, ...
    f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);



%% sort data by eKurtosis to help with the representation.

[dataForICA, sortedKurtDataIn, idxKurt] = sortDataByKurtosis(dataForICA);

signalNamesICA     = lisftOfSignals(idxKurt);
numInputSignalsICA = numel(signalNamesICA);

% Save some input info at this point
currentResults.point      = currentPoint;
currentResults.inputNames = signalNamesICA;
currentResults.eKurtInput = sortedKurtDataIn;
   
if doRemoveBaselineWandering
    fprintf('Removing baseline wandering...');
    for idx=1:numInputSignalsICA
        dataForICA(idx, :) = removeBaselineWandering(dataForICA(idx, :));
    end
end


%% Perform Botterom-Smith analysis.

% Get Botterom's frequency
[inputHatPSD, inputHatPSDOmega] = ...
    getPowerSpectralDensityBotteron(dataForICA, dataIn.fs, ...
    'windowType'   , 'rectangular', ...
    'windowSize'   , lengthSignal/4, ...
    'lengthFFT'    , 8192, ...
    'numberOverlap', [], ...
    'psdCutOff'    , 60);
 

% force to be between 4 and 9 Hz
idxValidFrequencies = inputHatPSDOmega >= 4;
idxValidFrequencies = idxValidFrequencies & inputHatPSDOmega < 10;

inputHatPSD(~idxValidFrequencies) = 0;


[inputHatFreqPeak, inputHatFreqPeakIdx] = max(inputHatPSD);
maxFrequency = inputHatPSDOmega(inputHatFreqPeakIdx);
botteronFrequency(idx) = maxFrequency;


% dumpdata
pointResults = [];
pointResults.pointId             = currentPoint;
pointResults.dataForICA          = dataForICA;
pointResults.inputHatPSD         = inputHatPSD;
pointResults.inputHatPSDOmega    = inputHatPSDOmega;
pointResults.inputHatFreqPeak    = maxFrequency;
pointResults.inputHatFreqPeakIdx = inputHatFreqPeakIdx;
cartoPointAnalysis{1, idx}       = pointResults;

end

srcDataName    = regexprep(srcDescription, {'\\', '\/'}, '_') ;
currentMatFile = ['Results', srcDataName, '.mat'];
vars2save      = {'cartoPointAnalysis', 'botteronFrequency', ...
                  'meshData', 'cartoData'};
save(currentMatFile, vars2save{:});

%%
figure,
xAxis = 1:max(allPoints);
botteronResults = zeros(size(xAxis));
botteronResults(allPoints) = botteronFrequency;
stem(xAxis, botteronResults);
title('Frequency Profile');
ylabel('Frequency (Hz)');
xlabel('Point Id');
%%
%cartoAnalysis = load(currentMatFile)
plotAnalysisOverMeshInteractive( meshData, cartoData, cartoPointAnalysis, botteronFrequency );

