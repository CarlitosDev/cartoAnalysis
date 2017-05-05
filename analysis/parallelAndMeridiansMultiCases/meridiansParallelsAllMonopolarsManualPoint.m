%% Set parameters
    testType   = 'plain (no filters, no baseline removal)';

    % set the signals to analyse
    pentaArray = listOfPentaArraySignals();

    referenceSignalName = pentaArray.reference{1};
    
    % Pass the set of unipolars and the algorithm will plot them
    signalsToResearch   = pentaArray.pentaUnipolar;
    
% Where to save results
    resultsFolder = fullfile(rootFolder(), 'resultsHighDensity');

% Path to CARTO file
    cartoFile = fullfile(getDataFolder(), 'RyC_Patient 2016_12_20', 'ECG.mat');

% Enable preprocess: filtering and baseline removal
    doPreprocess = false;
    doRemoveBaselineWandering = false;

%   Set the sampling frequency    
    samplingFrequency = 1e3;

% save processed file for visualisation
    doSaveProcessed = false;
    %doSaveProcessed = true;
    
% this is the point to analyse. 
    selectedPoint = 22;
    
% signal to analyse
  targetSignal = '20A_1';
  targetElectrode = 1;