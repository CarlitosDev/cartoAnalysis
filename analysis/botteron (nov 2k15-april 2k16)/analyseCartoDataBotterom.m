%% analyseCartoDataBotterom
% 
% This script runs an independent component analysis -fastICA- on ECG/EGM
% data given a LISFTOFSIGNALS and a LISTOFPOINTS.
% 
% To set up the script, follow the section below.
% 
%
%%

%listOfSelectedIcs = {[1,2]};
listOfSelectedIcs = {1};
%listOfSelectedIcs = {1, 2, 3, [1,2], [1,2,3]};

testType = 'NewDataS2m3m4';

% carto file (also relative to main folder)
%     path2CartoFile     = fullfile('data','cartoXP','handles.mat');   
%     fullPath2CartoFile = fullfile(mainFolder, path2CartoFile);
%path2CartoFile = fullfile('data','Auriculas','#4', 'Carto2', 'EG.mat');
%/Users/carlosAguilar/Documents/Matlab PhD/codigoCarlos/data/Auriculas/#5/Carto/1-Map
%path2CartoFile = fullfile('data','Auriculas','#5', 'Carto', '1-Map', 'EG.mat');
%path2CartoFile = fullfile('data','Auriculas','#4', 'Carto1', 'EG.mat');
path2CartoFile = fullfile('data','Auriculas','#2', 'Carto', 'EG.mat');
fullPath2CartoFile = fullfile(rootFolder(), path2CartoFile); 


numFlyBys = 1;
doRemoveBaselineWandering = false;

% Do whitening and filtering if requested    
    doPreprocess = true;

    % define input signals
    lisftOfSignals = {'M1_MINUS_M2', 'M3_MINUS_M4'};
    
    
    superiores = [113   107     8   137   108   115   383   124     5   308];
    inferiores = [287   272    62   264    99   286    40    26    88    33];
    medios     = [365   261   363   362   249   358   234   359    39   165];
    %listOfPoints = [superiores, inferiores, medios];
    %listOfPoints = 112:140;
    listOfPoints = 113;
    %listOfPoints = 'all';
        
    
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
    fid = fopen(fileName , 'a');
    
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
testResults    = cell(numSelectedICs, numPoints);


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

selDataIn = dataIn.signalData(idxSignals, :);

% Preprocess
% Schilling (p145) sets filters to (30, 150)
[numSignals, lengthSignal] = size(selDataIn);

f1HighPassCutOff = 1.0;
f2StopBandCutOff = 50 ;
f3LowPassCutOff  = 250;
dataForICA = preprocessForICAv2(dataIn, ...
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


fprintf('\nInput signals for fastICA:\n')
for idx=1:numInputSignalsICA
    fprintf('%d - %s\n', idx, signalNamesICA{idx});
end







end