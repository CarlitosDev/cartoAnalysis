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
    
% Do whitening and filtering if requested    
    doPreprocess = true;

    % define input signals
    lisftOfSignals = {'R1','R2','M1','M2','AVF','I','V1','V4','V6'};
    
    % set list of points here or 'all' for all points analysis.
    superiores = [113   107     8   137   108   115   383   124     5   308];
    inferiores = [287   272    62   264    99   286    40    26    88    33];
    medios     = [365   261   363   362   249   358   234   359    39   165];
    listOfPoints = [superiores, inferiores, medios];
    listOfPoints = 26;
    
    

%   Add paths for root, results and utils folders. Also set fastICA path
%   and the mat file where the data is stored.

    % main folder (work, vaio, mac...) 
    %mainFolder    = '/Users/carlosAguilar/Documents/Matlab PhD/codigoCarlos';
%    mainFolder    = 'C:\Colours\Matlab PhD\codigoCarlos';
    mainFolder    = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos';
    
    
    resultsFolder = fullfile(mainFolder, 'tests', 'signalPreProc');
    %mkdir(resultsFolder);

    % add utils to path
    currentFolder = 'utils';
    addpath(fullfile(mainFolder,currentFolder));

    % add fastICA     
    fastICAFolder = fullfile('external code', 'FastICA_2.5');
    addpath(fullfile(mainFolder, fastICAFolder));

    % move to tests folder
    currentFolder = 'tests';
    cd(fullfile(mainFolder, currentFolder));

    % carto file (also relative to main folder)
    path2CartoFile = fullfile('data','cartoXP','handles.mat');
    
    

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
    if listOfPoints(~presenPoints)
        fprintf('Warning: Point %d not found\n', listOfPoints(~presenPoints));
    end
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
if doPreprocess
    dataOut = preprocessForICA (dataIn);
else
  dataOut.signalData = dataIn.signalData;
end


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
addOne = @(x) x+1;
addsubplotindex = ...
  @(f,n) setappdata(f, 'subplotindices', [getappdata(f,'subplotindices'), n]);

numPlotColumns = 7;
numPlots       = 5;
offsetCol      = [1,0,1,0,0];
handlesArray   = zeros(1, numPlots*numOutputSignalsICA);

currentIdx     = 1;
currentHandle  = 1;
currentFig     = figure();

for idx = 1:numOutputSignalsICA  
  for idxCol = 1:numPlots      
  handlesArray(currentHandle) = subplot( numOutputSignalsICA, numPlotColumns, currentIdx + [0,offsetCol(idxCol)] );
  currentIdx    = 1 + currentIdx + offsetCol(idxCol);
  currentHandle = addOne(currentHandle);
  end
end

%% survey ICA output and plot

normalise = @(x) x./sum(x);

icaFeatures   = struct;
currentHandle = 1;

for idx=1:numOutputSignalsICA

  outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
  
  h0 = handlesArray(currentHandle); 
  plot (h0, dataForICA(idx,:), 'LineWidth', 2);
  title(h0, signalNamesICA{idx});
  axis (h0,'tight');

  currentHandle = addOne(currentHandle);
  h1 = handlesArray(currentHandle);
  
  % get the PSD for the current input signal
  [inputPSD, inputPSDOmega] = getPowerSpectralDensity(dataForICA(idx,:), dataIn.fs, ...
  'windowType', 'rectangular', 'windowSize', length(dataForICA(idx,:)), ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 30);

  inputPSD = normalise(inputPSD);  

  plot(h1, inputPSDOmega, inputPSD, 'LineWidth', 2);
  [~, inputFreqPeakIdx] = max(inputPSD);
  axis(h1, 'tight');
  hold(h1);
  stem(h1, inputPSDOmega(inputFreqPeakIdx), inputPSD(inputFreqPeakIdx), ...
  'LineWidth', 2);
  psdText  = sprintf('PSD F_d = %3.2f Hz', inputPSDOmega(inputFreqPeakIdx));
  title(h1, psdText);   
  set(h1, 'xtick',[]);
  set(h1, 'xticklabel',[]);
  set(h1, 'ytick',[]);
  set(h1, 'yticklabel',[]);   
  
 
  currentHandle = addOne(currentHandle);
  h2 = handlesArray(currentHandle);
  
  currentHandle = addOne(currentHandle);
  h3 = handlesArray(currentHandle);
  currentHandle = addOne(currentHandle);
  h4 = handlesArray(currentHandle); 
  
  [icaFeatures(idx).inputDataStats, icaFeatures(idx).powSpectDens, ...
    icaFeatures(idx).psdOmega, icaFeatures(idx).psdPeaks] = ...
    surveyIndependentComponents( icaSorted(idx, :), outSignalName, dataIn.fs, h2, h3, h4);
  
  set(h2,'ButtonDownFcn',@(~,~) addSubplotIndex(currentFig, idx, h2));
  currentHandle = addOne(currentHandle);
    
end

h = uicontrol('String', 'ICs selected', 'Position', [20 20 100 30], ...
'Callback', 'set(gcbf, ''Name'', sprintf(''IC %d '', getappdata(gcbf,''subplotindices'')))');
% this one doesn't get the figure 'completely' maximised
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% set(get(handle(currentFig),'JavaFrame'), 'Maximized', 1);
waitfor(currentFig, 'Name');

%end % remove that one

%% Re-construct signals with selected components
%
%  In this step we set to zero all of the sources regarded as VA and
%  reconstruct observations as $hat{x} = A \times s$

selectedICs = getappdata(currentFig, 'subplotindices');

% Get the observations from the sources. 
%   Set to zero all the ICs, and add back the ones picked.
%   Keep the original order as they were displatyed by eKurt order.
selectedSources = zeros(size(signalICA));
selectedSources(idxKurt(selectedICs), :) = icaSorted(selectedICs, :);

%   $hat{x} = A \times s$
dataForICAHat =  A * selectedSources;


% Overlay Xhat obsevations
currentHandle = 1;

for idx=1:numOutputSignalsICA
    
  h0 = handlesArray(currentHandle); 
  hold(h0); %hold current data
  plot(h0, dataForICAHat(idx,:), 'LineWidth', 1, 'Color', 'r');
  
  currentHandle = addOne(currentHandle);
  h1 = handlesArray(currentHandle); 
  [inputHatPSD, inputHatPSDOmega] = getPowerSpectralDensity(dataForICAHat(idx,:), dataIn.fs, ...
  'windowType', 'rectangular', 'windowSize', length(dataForICAHat(idx,:)), ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 30);
  
  % normalise the PDS so we can overlay the Xhat PSD
  inputHatPSD = normalise(inputHatPSD);  

  hold(h1, 'on'); %hold current data
  plot(h1, inputHatPSDOmega, inputHatPSD, 'LineWidth', 1, 'Color', 'r');
  [~, inputHatFreqPeakIdx] = max(inputHatPSD);
  hold(h1, 'on');
  stem(h1, inputHatPSDOmega(inputHatFreqPeakIdx), inputHatPSD(inputHatFreqPeakIdx), ...
  'LineWidth', 1, 'Color', 'r');
  axis(h1, 'tight');

  currentPSDText = get(get(h1, 'title'),'String');% get text back before overwritting
  psdText = sprintf('%s\nIC F_d = %3.2f Hz', currentPSDText, inputHatPSDOmega(inputHatFreqPeakIdx));
  title(h1, psdText);
  
  currentHandle = currentHandle + numPlots-1;  
  
end

% update results
for idx = selectedICs
    icaFeatures(idx).icSelected = 1;
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