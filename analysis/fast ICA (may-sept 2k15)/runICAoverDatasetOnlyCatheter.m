% %% Run fastICA over a collection of signals
% 
% This script runs an independent component analysis -fastICA- on ECG/EGM
% data given a LISFTOFSIGNALS and a LISTOFPOINTS.
% 
% To set up the script, follow the section below.
% 
% The OUTPUT of the script -results folder- will be:
%
%     - A collection of png files named <Point XX.png> where XX stands for
%     the current point. The images summarise the input signals, the
%     independent components and their PSD and histogram. 
%
%     - A mat file <icaTestDDDD.mat> where DDDD is the creation date. The
%     contents of the file as organised as follows:
%
%         * The analysis of each point is stored in a cell and can be
%         directly accessed by its index.
%         * Every point has a set of structs that contain the measures of
%         the IC - sorted by kurtosis in ascend order.
%         * The number of structs depends on the number of IC's that
%         fastICA managed to produce.
%         * Every struct contains statistical measures in INPUTDATASTATS
%         and PSD-derived measures: PSD itself, PSD frequency axis and the
%         main peaks of the PSD.
% 
% 
%       Reading the output. Examples:
%
%       (1) To get the fd of the IC with lower kurtosis of the point 84 :  
%           fd = testResults{1,84}(1,1).psdPeaks(1);
% 
%       (2) To get the kurtosis of all ICs for the 48th:
%           for idx=1:numel(testResults{1,48})
%             fprintf(' IC %d - eKurtosis %3.2f\n', idx, ...
%                testResults{1,48}(1,idx).inputDataStats.eKurtosis);
%           end
% 
% 
%
%                                                  SVN fields:
%
%                                                  $Revision: $
%                                                  $Author: carlosAguilar $
%                                                  $Date: $

%% Set up configuration

%   Set the sampling frequency, the list of input signals for fastICA and
%   the points to extract the data from. The latter can be just a list of
%   points or the string 'all' meaning that the process will be repeated
%   through all the data.

    variables.fs = 1000;

    % define input signals
    lisftOfSignals = {'R1','R2','M1','M2'};
    
    % set list of points here or 'all' for all points analysis.
    listOfPoints = [1,48,84,100];
    
    

%   Add paths for root, results and utils folders. Also set fastICA path
%   and the mat file where the data is stored.

    % main folder
    mainFolder    = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos';   
    
    
    resultsFolder = fullfile(mainFolder, 'tests\allPoints');

    % add utils to path
    currentFolder = 'utils';
    addpath(fullfile(mainFolder,currentFolder));

    % add fastICA     
    fastICAFolder = 'external code\FastICA_2.5';
    addpath(fullfile(mainFolder, fastICAFolder));

    % move to tests folder
    currentFolder = 'tests';
    cd(fullfile(mainFolder, currentFolder));

    % carto file (also relative to main folder)
    path2CartoFile = 'data\cartoXP\handles.mat';
    
    

%% No need to edit below this point
% -----------------------------------


% load CartoXP
load(fullfile(mainFolder, path2CartoFile));
information  = handles.DB;
data         = information.EG;


% Set RNG for reproducibility, keeping the old settings to be
% restored at the end. 
oldRNG = RandStream.getGlobalStream;
cleaners.RNG =...
  onCleanup(@()RandStream.setGlobalStream(oldRNG));

newRNG = RandStream('mrg32k3a', 'Seed', 0);
RandStream.setGlobalStream(newRNG);

%% where to apply fastICA.

if isa(listOfPoints, 'char') && strcmpi(listOfPoints, 'all')
    allPoints = unique(information.npunto);  
    % let's assume numeric otherwise    
else
    % check the requested point belongs to the list. Otherwise, just wipe
    % it out.  
    currentPoints = unique(information.npunto);
    presenPoints  = ismember(listOfPoints, currentPoints);
    allPoints     = listOfPoints(presenPoints);   
    fprintf('Warning: Point %d not found\n', listOfPoints(~presenPoints));
end    

numPoints    = numel(allPoints);
testResults  = cell(1, numPoints);


for idxPoint = 1:numPoints
  
currentPoint = allPoints(idxPoint);
  
idx      = find(information.npunto == currentPoint);
ecgData  = data(idx);
figTitle =  sprintf('Point %d', currentPoint);

%% Preprocess input data. Be verbose
dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;

for idx=1:numel(dataIn.signalNames)
    fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
end

dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;
dataOut            = preprocessForICA (dataIn);


%% Check all signals are present

numInputSignalsICA  = numel(lisftOfSignals);
[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');
dataForICA     = dataOut.signalData(idxSignals, :);
signalNamesICA = dataIn.signalNames(idxSignals);


%% run fastICA and sort output by excess Kurtosis
[signalICA, A, W] = fastica(dataForICA);

kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);
icaSorted = signalICA(idxKurt, :);

numOutputSignalsICA = numel(kurtosisICA);

if numInputSignalsICA ~= numOutputSignalsICA
  warning('fastICA found %d independent components (%d signals slashed)\n', ...
  numOutputSignalsICA, numInputSignalsICA- numOutputSignalsICA); %#ok<WNTAG>
end


%% prepare data for representation
numPlotColumns = 6;
numPlots       = 4;
offsetCol      = [1,1,0,0];
handlesArray   = zeros(1, numPlots*numOutputSignalsICA);

currentIdx     = 1;
currentHandle  = 1;
currentFig     = figure();

for idx = 1:numOutputSignalsICA  
  for idxCol = 1:4       
  handlesArray(currentHandle) = subplot( numOutputSignalsICA, numPlotColumns, currentIdx + [0,offsetCol(idxCol)] );
  currentIdx = 1 + currentIdx + offsetCol(idxCol);
  currentHandle = currentHandle + 1;
  end
end

%% survey ICA output and plot
addOne = @(x) x+1;

icaFeatures   = struct;
currentHandle = 1;

for idx=1:numOutputSignalsICA

  outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
  
  h0 = handlesArray(currentHandle); 
  plot (h0, dataForICA(idx,:))
  title(h0, signalNamesICA{idx});
  axis (h0,'tight')
  
  currentHandle = addOne(currentHandle);
  h1 = handlesArray(currentHandle);
  currentHandle = addOne(currentHandle);
  h2 = handlesArray(currentHandle);
  currentHandle = addOne(currentHandle);
  h3 = handlesArray(currentHandle); 
  
  [icaFeatures(idx).inputDataStats, icaFeatures(idx).powSpectDens, ...
    icaFeatures(idx).psdOmega, icaFeatures(idx).psdPeaks] = ...
    surveyIndependentComponents( icaSorted(idx, :), outSignalName, dataIn.fs, h1, h2, h3);
  
  currentHandle = addOne(currentHandle);
    
end

% save figure
set(currentFig, 'units', 'normalized', 'outerposition',[0 0 1.5 1.5]);
set(currentFig, 'Name' , figTitle);

figPathName = fullfile(resultsFolder, [figTitle, '.png']);
hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
close(currentFig);
fprintf('\n');

% 
testResults{currentPoint} = icaFeatures;

end

save(fullfile(resultsFolder, ['icaTest', datestr(now, 'dd-mm-yyyy hh'), 'h.mat']), 'testResults');