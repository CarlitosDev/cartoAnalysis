%% testMeridiansAndParallels
%
% This script loads the cartoXP data stored in .\Data (ungitted)\Auriculas
% using 'loadCartoXPData'
% 
%
%% Set parameters

testType   = 'simBotteronSmith';

% Where to save results
    destFolder = fullfile(rootFolder(), 'results');


% CARTO file
    cartoFile = fullfile('#2', 'Carto', 'EG.mat');
    %cartoFile = fullfile('#4', 'Carto1', 'EG.mat');
    %cartoFile = fullfile('#5', 'Carto', '1-Map', 'EG.mat');    

% Enable preprocess: filtering and baseline removal
    doPreprocess = true;
    doRemoveBaselineWandering = true;

% save processed file for visualisation
  doSaveProcessed = true;

%% Set up configuration

%   Set the sampling frequency, the list of input signals for fastICA and
%   the points to extract the data from. The latter can be just a list of
%   points or the string 'all' meaning that the process will be repeated
%   through all the data.

    variables.fs = 1000;   
    
    inputSignalName     = 'M1_MINUS_M2';
    referenceSignalName = 'V1';
    listOfSignals       = {inputSignalName, referenceSignalName};   

    superiores = [113   107     8   137   108   115   383   124     5   308];
    inferiores = [287   272    62   264    99   286    40    26    88    33];
    medios     = [365   261   363   362   249   358   234   359    39   165];
    %listOfPoints = [superiores, inferiores, medios];
    %listOfPoints = 112:140;
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


%% load CartoXP (same interface for loadSimulation)
[data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile);

%% Copy the analysis points selected above

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


    fprintf('\nRunning analysis on signal %d\n', currentPoint);

    %% Check all signals are present and then preprocess input data.
    dataIn = [];
    dataIn.signalNames = ecgData.tipo_ECG;
    dataIn.signalData  = cell2mat(ecgData.signal');
    dataIn.fs          = variables.fs;

    [signalPresent, idxSignals] = ismember(listOfSignals, dataIn.signalNames);
    assert(all(signalPresent), 'Missing signals');

    selDataIn.signalData  = dataIn.signalData(idxSignals, :);
    selDataIn.signalNames = dataIn.signalNames(idxSignals);
    selDataIn.fs          = dataIn.fs;

    %% Filter preprocess
    if doPreprocess
        % Schilling (p145) sets filters to (30, 150)
        [numSignals, lengthSignal] = size(selDataIn.signalData);
        dataForICA = preprocessForICAv2(selDataIn, ...
            f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);
    else
        dataForICA  = selDataIn.signalData;
    end

%% Once the signals are filtered, split into reference and analysis

    [~, refSignalIdx] = ismember(referenceSignalName, selDataIn.signalNames);
    referenceSignal = dataForICA(refSignalIdx, :);
    % flush it
    dataForICA(refSignalIdx, :) = [];
    
    %plotEGMandECGTimeAndSpectrum(dataForICA(1, :), referenceSignal, 1e3);

%% sort data by eKurtosis to help with the representation.

    [dataForICA, sortedKurtDataIn, idxKurt] = sortDataByKurtosis(dataForICA);

    signalNamesICA     = listOfSignals(idxKurt);
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
%     getPowerSpectralDensityBotteron(dataForICA, dataIn.fs, ...
%     'windowType'   , 'rectangular' , ...
%     'windowSize'   , lengthSignal/4, ...
%     'lengthFFT'    , 8192, ...
%     'numberOverlap', [], ...
%     'psdCutOff'    , 60);


% get the PSD for the input signal
[inputPSD, inputPSDOmega] = getPowerSpectralDensity(dataForICA, dataIn.fs, ...
'windowType', 'rectangular', 'windowSize', length(dataForICA), ...
'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 60);
 
%% Use JL's to get main frequencies. Name this signal 'hat'

[maxFrequency,indff,bwfromf0,z,periodogramInput,fz] = df_Ng(dataForICA, dataIn.fs);
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

% Cut-off freqs above 35 Hz
idxValidFrequencies = fz <= 35;

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

end

%%

if doSaveProcessed 
  currentMatFile = ['Results ', srcDataName, '.mat'];
  vars2save      = {'cartoPointAnalysis', 'botteronFrequency', ...
                    'meshData', 'cartoData'};
  save(currentMatFile, vars2save{:});
end



%% Check the new visualisation
%plotMeshAndAllEGMAnalysis( meshData, cartoData, cartoPointAnalysis, botteronFrequency)
plotMeshAndAllEGMAnalysisParallel( meshData, cartoData, cartoPointAnalysis, botteronFrequency)