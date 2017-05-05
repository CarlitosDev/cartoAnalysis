%% meridiansParallelsPredefinedPoint
%
% This script loads the new cartoXP data stored in .\Data (ungitted)\FA_RamonYCajal
% using 'loadRamonyCajalData'
% Instead of showing a figure that allows to select a point, the point
% number is hard-coded below. This method is useful to produce the analysis plots
% using different signals.
%
%   Carlos Aguilar - 16/02/2017
%

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
    
   

%% Set up configuration

    listOfPoints = 'all';
    
% ica parameters
    numFlyBys = 1;
    listOfSelectedIcs = {1};
    % to group outputs:
    %listOfSelectedIcs = {1, 2, 3, [1,2], [1,2,3]};
    
    
% Preprocessing features
    f1HighPassCutOff = 1.0;
    f2StopBandCutOff = 50 ; % notch
    f3LowPassCutOff  = 250;


%% work path out as it is ungitted

    matchingFNames = which(cartoFile);
    if isempty(matchingFNames) 
        errordlg('Can''t resolve path to file. Run initialise()');
        return;
    end
    fullPath2CartoFile = matchingFNames;
    srcDataName  = regexprep(cartoFile, {'\\', '\/'}, ' ');
        
    
%%   Add paths for root, results and utils folders. Also set fastICA path
%   and the mat file where the data is stored.

    commaString = @(listNumbers) ...
      strcat(sprintf('%d,', listNumbers(1:end-1)), sprintf('%d', listNumbers(end)));   
    signalsStr  = @(listSignals) ...
      [sprintf('%s,', listSignals{1:end-1}), sprintf('%s', listSignals{end})];
    floatStr    = @(listNumbers) ...
      strcat(sprintf('%3.3f,', listNumbers(1:end-1)), sprintf('%3.3f', listNumbers(end)));

    fileName = fullfile(resultsFolder, ['log', datestr(now), '.csv']);
    fid = 1;
    %fid = fopen(fileName , 'a');


%% load loadRamonyCajalData (same interface for loadSimulation/loadCartoXPData)

[data, pointsInfo, meshData, cartoData, cathererData] = ...
  loadRamonyCajalData(fullPath2CartoFile);

cathTable    = struct2table(cathererData);
idxElectrode = cathTable.electrodeID == targetElectrode;
cathTable    = cathTable(idxElectrode, :);
cathTable.pointsId(:) = 1:height(cathTable);

pointsPosition = cathTable{:, {'x', 'y', 'z'}};

cathererData = table2struct(cathTable);

%   processedDataFolder = fullfile(rootFolder(), '..', 'Data (ungitted)', ...
%                         'preprocessed data');
%   currentMatFile      = ['raw_', srcDataName, '.mat'];
%   fullPath2MatFile    = fullfile(processedDataFolder, currentMatFile);
% 
%   vars2save           = {'data', 'pointsInfo', ...
%                          'meshData', 'cartoData', 'cathererData'};
%   
%   save(fullPath2MatFile, vars2save{:});

%load('D:\carlosAguilar\PhD Carlos\Data (ungitted)\FA_RamonYCajal\rawData_FA_RamonYCajal_ECG.mat')

%   processedDataFolder = fullfile(rootFolder(), '..', 'Data (ungitted)', ...
%                         'FA_RamonYCajal');
%   currentMatFile      = ['rawData_FA_RamonYCajal_ECG', '.mat'];
%   fullPath2MatFile    = fullfile(processedDataFolder, currentMatFile);  
%   load(fullPath2MatFile);

% fprintf('\n Delete this bit that overrides the current listOfSignals\n')
% listOfSignals = data.tipo_ECG;

validPoint = pointsInfo.pointToIdxMap.isKey(selectedPoint);
assert(validPoint, 'Point %d does not exist', selectedPoint);

%% Get all the points from the catherers

    if isa(listOfPoints, 'char') && strcmpi(listOfPoints, 'all')
        allPoints = [cathererData.pointsId];
        % let's assume numeric otherwise    
    else
        % check the requested point belongs to the list. Otherwise, just wipe
        % it out.  
        currentPoints = pointsInfo.allPointsInfo;
        presenPoints  = ismember(listOfPoints, currentPoints);
        allPoints     = listOfPoints(presenPoints);
        if listOfPoints(~presenPoints)
            fprintf('Warning: Point %d not found\n', listOfPoints(~presenPoints));
        end
    end

    numPoints      = numel(allPoints);
    numSelectedICs = numel(listOfSelectedIcs);
    testResults    = cell(numSelectedICs, numPoints);


%     % Populate the header
%     fprintf(fid, '%s\t%s\t%s\t%s\ts%s\t%s\n', ...
%         'Point', 'Input signals', 'input eKurtosis', ...
%         'number of ICs', 'output eKurtosis', 'output peak Freq');


for idxPoint=1:numPoints
    
    currentResults = [];
    currentPoint   = allPoints(idxPoint);
    idx = cathererData(currentPoint).uniquePentaRayID;
    pentaRayElectrode = cathererData(currentPoint).pentaRayElectrode;
    
    fprintf('\nRunning Pentaray analysis on point %d\n', currentPoint);
    
    allFNames = fieldnames(data);
    % some signals are 'data' instead of 'signal'...
    aliasData = {'signal', 'data'};
    dataFName = aliasData{ismember(aliasData, allFNames)};
    
    
    %% Look for the reference signal
    
    idxReference = contains(data(idx).tipo_ECG, referenceSignalName);
    ecgRefData   = data(idx).(dataFName){idxReference};
    ecgGain      = data(idx).gain;

    %% Check all signals are present and then preprocess input data.
    
    selDataIn = [];
    selDataIn.signalNames = {referenceSignalName, pentaRayElectrode};
    signalData         = ...
      vertcat(ecgRefData, cathererData(currentPoint).signal);
    selDataIn.signalData  = signalData.*ecgGain;
    selDataIn.fs          = samplingFrequency;
    
    [numSignals, lengthSignal] = size(selDataIn.signalData);

    %% Filter preprocess
    if doPreprocess
        % Schilling (p145) sets filters to (30, 150)
        dataForICA = preprocessForICAv2(selDataIn, ...
            f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);
    else
        dataForICA  = selDataIn.signalData;
    end

%% Once the signals are filtered, split into reference and analysis

    [~, refSignalIdx] = ismember(referenceSignalName, selDataIn.signalNames);
    referenceSignal   = dataForICA(refSignalIdx, :);
    % flush it
    dataForICA(refSignalIdx, :) = [];
    
    %plotEGMandECGTimeAndSpectrum(dataForICA(1, :), referenceSignal, 1e3);

%% sort data by eKurtosis to help with the representation.
    [dataForICA, sortedKurtDataIn, idxKurt] = sortDataByKurtosis(dataForICA);

    signalNamesICA     = selDataIn.signalNames(idxKurt);
    numInputSignalsICA = numel(signalNamesICA);

    % Save some input info at this point
    currentResults.point      = currentPoint;
    currentResults.inputNames = signalNamesICA;
    currentResults.eKurtInput = sortedKurtDataIn;

    if doRemoveBaselineWandering
        % sample the original signal for the b-spline in chunks of 20% the total size
        samplingStep = 20;
        fprintf('Removing baseline wandering...');
        for idx=1:numInputSignalsICA
            [dataForICA(idx, :), ~] = ...
                removeBaselineWanderingV2(dataForICA(idx, :), samplingStep);        
        end
    end


%% Perform Botterom-Smith analysis.

% Get Botterom's frequency
% [inputHatPSD, inputHatPSDOmega] = ...
%     getPowerSpectralDensityBotteron(dataForICA, selDataIn.fs, ...
%     'windowType'   , 'rectangular' , ...
%     'windowSize'   , lengthSignal/4, ...
%     'lengthFFT'    , 8192, ...
%     'numberOverlap', [], ...
%     'psdCutOff'    , 60);


% get the PSD for the input signal
[inputPSD, inputPSDOmega] = getPowerSpectralDensity(dataForICA, selDataIn.fs, ...
'windowType', 'rectangular', 'windowSize', length(dataForICA), ...
'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 60);
 
%% Use JL's to get main frequencies. Name this signal 'hat'

[maxFrequency,indff,bwfromf0,z,periodogramInput,fz] = df_Ng(dataForICA, selDataIn.fs);
% deltaFreq = fz(2)-fz(1);
% ancho = 2*deltaFreq;
ancho = 1.8;%1Hz
[oi,ejexplot,ejeyplot] = organizationIndex(periodogramInput,fz,maxFrequency,ancho);
fprintf('Frequency Peak     %3.2f\n', maxFrequency);
fprintf('Organisation Index %3.2f\n', oi);

% % force to be between 4 and 9 Hz
% idxValidFrequencies = inputHatPSDOmega >= 4;
% idxValidFrequencies = idxValidFrequencies & inputHatPSDOmega < 10;
% 
% inputHatPSD(~idxValidFrequencies) = 0;

% Cut-off freqs above 100 Hz
idxValidFrequencies = fz <= 100;

inputHatFreqPeakIdx = indff;
inputHatPSD = periodogramInput(idxValidFrequencies);
inputHatPSDOmega = fz(idxValidFrequencies);
botteronFrequency(idxPoint) = maxFrequency;


% dumpdata
pointResults = [];
pointResults.pointId             = currentPoint;
% input signal
pointResults.dataForICA          = dataForICA;
pointResults.inputPSD            = inputPSD;
pointResults.inputPSDOmega       = inputPSDOmega;
% analysed signal ('hat')
pointResults.inputHatPSD         = inputHatPSD;
pointResults.inputHatPSDOmega    = inputHatPSDOmega;
pointResults.inputHatFreqPeak    = maxFrequency;
pointResults.inputHatFreqPeakIdx = inputHatFreqPeakIdx;
pointResults.organisationIndex   = oi;
pointResults.searchWindow        = ancho;
pointResults.referenceSignal     = referenceSignal;
cartoPointAnalysis{1, idxPoint}  = pointResults;


% plotSignalAndSpectrum(dataForICA, selDataIn.fs);
% drawnow
% 
%  figTitle = sprintf('%s_point_%d', ...
%             inputSignalName, currentPoint);
%  title(figTitle);
%  saveFigureAsPNG(figTitle, 'D:\carlosAguilar\PhD Carlos\tempPlots');
% 
%  close(gcf());

end

%%

if doSaveProcessed 
  processedDataFolder = fullfile(rootFolder(), '..', 'Data (ungitted)', ...
                        'preprocessed data');
  currentMatFile      = ['PreProcessed_', srcDataName, '.mat'];
  fullPath2MatFile    = fullfile(processedDataFolder, currentMatFile);

  vars2save           = {'cartoPointAnalysis', 'botteronFrequency', ...
                         'meshData', 'cartoData'};

  save(fullPath2MatFile, vars2save{:});
end


%% Prepare the data for the autovisualisation

    % Let's assume that the atrium is quite round. Get the centroid and
    % the direction of a vector from the centroid to every point. Then,
    % shift the text by 'delta' pointing outwards.
    deltaText      = 3.0;
    
    pointsCentroid = mean(pointsPosition, 1);
    pointsCentered = ...
        bsxfun(@minus, pointsPosition, pointsCentroid);
    
    vectorDirection = sign(pointsCentered);
    
    textPosition = pointsPosition + deltaText.*vectorDirection;
    
    % arrange carto data
    xData = pointsPosition(:,1);
    yData = pointsPosition(:,2);
    zData = pointsPosition(:,3);

    % Center points
    pointsCentroid = mean(pointsPosition, 1);
    pointsCentered = ...
        bsxfun(@minus, pointsPosition, pointsCentroid);

    % Move along the azimuth. Go to spherical coordinates to sort the points by
    % azimuth, so we can move the camera along the horizontal axis.
    [currentAzimuth, currentElevation, currentRadius] = ...
        cart2sph(pointsCentered(:,1), ...
        pointsCentered(:,2), ...
        pointsCentered(:,3));

    fprintf('\nSelected point for analysis is %d\n', selectedPoint);      
    selectedPointIdx = selectedPoint;
      
    [parallelPath, meridianPath, pathInfo] = ...
        getMeridianAndParallelV4(currentAzimuth, currentElevation, ...
        currentRadius, selectedPointIdx, pointsCentered);

    % Arrange selected data
    selectedPointsData = [];
    
    % add some info
      selectedPointsData.testType        = testType;
      selectedPointsData.inputSignalName = cartoFile;
      selectedPointsData.fs              = samplingFrequency;
      selectedPointsData.lengthSignal    = lengthSignal;
      
    % push the data from the analysis
      selectedPointsData.selectedPoint    = selectedPoint;
      selectedPointsData.selectedPointIdx = selectedPointIdx;
      selectedPointsData.selectedAz       = currentAzimuth(selectedPointIdx);
      selectedPointsData.selectedEl       = currentElevation(selectedPointIdx);
      selectedPointsData.textPosition     = textPosition;
      selectedPointsData.pointsId         = [cathererData.pointsId];      
      %selectedPointsData.pointsId         = cartoData.pointsId;

      
    % add summaries of the paths
      selectedPointsData.numParallelPoints    = pathInfo.numParallelPoints;
      selectedPointsData.allParallelPointsIdx = pathInfo.allParallelPointsIdx;      
      selectedPointsData.numMeridianPoints    = pathInfo.numMeridianPoints;
      selectedPointsData.allMeridianPointsIdx = pathInfo.allMeridianPointsIdx;

    % add the paths
      selectedPointsData.parallelPath         = parallelPath;
      selectedPointsData.meridianPath         = meridianPath;
      selectedPointsData.pathInfo             = pathInfo;

plotDataAnalysisParallelAndMeridianV2(cartoPointAnalysis, meshData, ...
  xData, yData , zData, selectedPointsData);