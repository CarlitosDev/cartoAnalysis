% %% Run fastICA over a collection of signals
% 
% This script runs an independent component analysis -fastICA- on ECG/EGM
% data given a LISFTOFSIGNALS and a LISTOFPOINTS.
% 
% To set up the script, follow the section below.
% 
%
%%

%selectedICs = [1,3];
selectedICs = [1,2,3];

% Do whitening and filtering if requested    
    doPreprocess = true;

    % define input signals
    lisftOfSignals = {'M1', 'M2', 'M3', 'M4', 'R1', 'R2', 'V1', 'V2', 'II'};
    
    % set list of points here or 'all' for all points analysis.
    superiores = [113   107     8   137   108   115   383   124     5   308];
    inferiores = [287   272    62   264    99   286    40    26    88    33];
    medios     = [365   261   363   362   249   358   234   359    39   165];
    listOfPoints = [superiores, inferiores, medios];
        


%%   Add paths for root, results and utils folders. Also set fastICA path
%   and the mat file where the data is stored.

    % main folder (work, vaio, mac...) 
    %mainFolder    = '/Users/carlosAguilar/Documents/Matlab PhD/codigoCarlos';
%    mainFolder    = 'C:\Colours\Matlab PhD\codigoCarlos';

    commaString = @(listNumbers) strcat(sprintf('%d,', listNumbers(1:end-1)), sprintf('%d', listNumbers(end)));

    setUp = false;

    mainFolder    = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos';
    
    % for reconstructing from ICA
    
    
    resultsFolder = fullfile(mainFolder, 'tests', ['ics_', commaString(selectedICs)]);
    mkdir(resultsFolder);

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


    
    

        
    
%% Set up configuration

%   Set the sampling frequency, the list of input signals for fastICA and
%   the points to extract the data from. The latter can be just a list of
%   points or the string 'all' meaning that the process will be repeated
%   through all the data.

    variables.fs = 1000;   
        

%% load CartoXP
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


figTitle =  sprintf('Point %d_ICs %s', currentPoint, commaString(selectedICs));

%% Preprocess input data. Be verbose
dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;

% for idx=1:numel(dataIn.signalNames)
%     fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
% end

dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;


%% Check all signals are present


[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

assert(all(signalPresent), 'Missing signals');

selDataIn = dataIn.signalData(idxSignals, :);

% remove the mean
[numSignals, lengthSignal]  = size(selDataIn);
selDataIn = selDataIn - mean(selDataIn,2)*ones(1,lengthSignal);


%% sort data by eKurtosis before ICA

kurtosisDataIn = kurtosis(selDataIn, 1, 2) - 3;
[~, idxKurt]   = sort(kurtosisDataIn);
selDataIn      = selDataIn(idxKurt, :);

signalNamesICA     = lisftOfSignals(idxKurt);
numInputSignalsICA = numel(signalNamesICA);


%% Preprocess here

%   TEST 1:
%
% The preprocessing consists of bandpass filtering (finite impulse response
% (FIR), 40–250 Hz, order 40, Kaiser window), rectification, and lowpass
% filtering (FIR, 0–20 Hz, order 40, Kaiser
% window)~cite{AFPreproc:richter2011novel}.
filterOrder    = 128;
dataSignalProc = filterFIR(selDataIn, dataIn.fs, ...
        'filterOrder'     , filterOrder, ...
        'filterType'      ,'high'      , ...
        'filterFrequency' , 1.0       , ...
        'showFilterDesign', false);
    
dataForICA = filterFIR(dataSignalProc, dataIn.fs, ...
        'filterOrder'     , filterOrder, ...
        'filterType'      ,'low'       , ...
        'filterFrequency' , 40         , ...        
        'showFilterDesign', false);
    


fprintf('\nInput signals for fastICA:\n')
for idx=1:numInputSignalsICA
    fprintf('%d - %s\n', idx, signalNamesICA{idx});
end


%% run fastICA and sort output by excess Kurtosis

%covariance = cov(dataForICA');

[signalICA, A, W] = fastica(dataForICA);

kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);
icaSorted = signalICA(idxKurt, :);

numOutputSignalsICA = numel(kurtosisICA);

if numInputSignalsICA ~= numOutputSignalsICA
  warning('fastICA found %d independent components (%d signals slashed)\n', ...
  numOutputSignalsICA, numInputSignalsICA- numOutputSignalsICA); %#ok<WNTAG>
end


%% Re-construct signals with selected components
%
%  In this step we set to zero all of the sources regarded as VA and
%  reconstruct observations as $hat{x} = A \times s$

% Get the observations from the sources. 
%   Set to zero all the ICs, and add back the ones picked.
%   Keep the original order as they were displatyed by eKurt order.
selectedSources = zeros(size(signalICA));
selectedSources(idxKurt(selectedICs), :) = icaSorted(selectedICs, :);

%   $hat{x} = A \times s$
dataForICAHat =  A * selectedSources;


%% prepare data for representation
addOne = @(x) x+1;

numPlotColumns = 6;
numPlots       = 4;
offsetCol      = [1,0,1,0];
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
  
  % plot time signals
      h0 = handlesArray(currentHandle); 
      plot(h0, selDataIn(idx, :), 'Color', [0 1 1], 'LineWidth', 2.5);  
      hold(h0);
      plot(h0, dataForICA(idx,:), 'LineWidth', 1.5);
      title(h0, signalNamesICA{idx});
      axis(h0,'tight');
      plot(h0, dataForICAHat(idx,:), 'LineWidth', 1.5, 'Color', 'r');
  
  
  currentHandle = addOne(currentHandle);
  h1 = handlesArray(currentHandle);
  
  % Plot PSDs
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
      psdText  = sprintf('F_d = %3.2f Hz', inputPSDOmega(inputFreqPeakIdx));
      title(h1, psdText);   
      set(h1, 'xtick',[]);
      set(h1, 'xticklabel',[]);
      set(h1, 'ytick',[]);
      set(h1, 'yticklabel',[]);

      [inputHatPSD, inputHatPSDOmega] = getPowerSpectralDensity(dataForICAHat(idx,:), dataIn.fs, ...
      'windowType', 'rectangular', 'windowSize', length(dataForICAHat(idx,:)), ...
      'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 30);

      inputHatPSD = normalise(inputHatPSD);

      plot(h1, inputHatPSDOmega, inputHatPSD, 'LineWidth', 2.0, 'Color', 'r');
      [~, inputHatFreqPeakIdx] = max(inputHatPSD);
      hold(h1, 'on');
      stem(h1, inputHatPSDOmega(inputHatFreqPeakIdx), inputHatPSD(inputHatFreqPeakIdx), ...
      'LineWidth', 1.5, 'Color', 'r');
      axis(h1, 'tight');
  
 
  currentHandle = addOne(currentHandle);
  h2 = handlesArray(currentHandle);
  
  currentHandle = addOne(currentHandle);
  h3 = handlesArray(currentHandle);

  [icaFeatures(idx).inputDataStats, icaFeatures(idx).powSpectDens, ...
    icaFeatures(idx).psdOmega, icaFeatures(idx).psdPeaks] = ...
    surveyIndependentComponents( icaSorted(idx, :), outSignalName, dataIn.fs, h2, h3, []);

  if any(ismember(idx, selectedICs))
      set(h2, 'LineWidth',2, 'Color', [1 1 0.9]);
  end
  
  currentHandle = addOne(currentHandle);
    
end

set(gcf,'units','normalized','outerposition',[0 0 1 1]);

% save figure
%set(currentFig, 'units', 'normalized', 'outerposition',[0 0 1.5 1.5]);
set(currentFig, 'Name' , figTitle);

figPathName = fullfile(resultsFolder, [figTitle, '.png']);
hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
close(currentFig);
fprintf('\n');
end