% %% Run fastICA over a collection of signals
% 
% This script runs an independent component analysis -fastICA- on ECG/EGM
% data given a LISFTOFSIGNALS and a LISTOFPOINTS.
% 
% To set up the script, follow the section below.
% 
%
%%

listOfSelectedIcs = {[1,2]};
%listOfSelectedIcs = {1, 2, 3, [1,2], [1,2,3]};

testType = 'Both';
numFlyBys = 1;
doRemoveBaselineWandering = true;

% Do whitening and filtering if requested    
    doPreprocess = true;

    % define input signals
    %lisftOfSignals = {'M1_MINUS_M2', 'M3_MINUS_M4', 'M1', 'M2', 'M3', 'M4', 'R1', 'R2', 'V1', 'V2', 'II'};
    %lisftOfSignals = {'M1_MINUS_M2', 'M3_MINUS_M4', 'M1', 'M2', 'M3', 'M4'};
    %lisftOfSignals = {'V1','V2','V3','V6', 'I','II'};
    %lisftOfSignals = {'M1', 'M2', 'M3', 'M4', 'R1', 'R2', 'V1', 'V2', 'II'};
    lisftOfSignals = {'M1', 'M2', 'M3', 'M4', 'R1', 'R2'};
    
    % set list of points here or 'all' for all points analysis.
    pointA = [107.8, 14.35, 2.257];
    pointB = [125.4, 53.9, -15.52];
    pointC = [133.7, 24.3, -50.48];

    
    
    
    superiores = [113   107     8   137   108   115   383   124     5   308];
    inferiores = [287   272    62   264    99   286    40    26    88    33];
    medios     = [365   261   363   362   249   358   234   359    39   165];
    %listOfPoints = [superiores, inferiores, medios];
    %listOfPoints = superiores
    listOfPoints = 'all';
        
    

 %% Set fastICA parameters in here
fastICA = [];
fastICA.approach = 'defl'; %'symm'; % 'defl'
fastICA.g = 'pow3'; %'tanh', 'gauss', 'skew','pow3'
fastICA.stabilization = 'off';

%%   Add paths for root, results and utils folders. Also set fastICA path
%   and the mat file where the data is stored.

    commaString = @(listNumbers) strcat(sprintf('%d,', listNumbers(1:end-1)), sprintf('%d', listNumbers(end)));   
    signalsStr  = @(listSignals) [sprintf('%s,', listSignals{1:end-1}), sprintf('%s', listSignals{end})];
    floatStr    = @(listNumbers) strcat(sprintf('%3.3f,', listNumbers(1:end-1)), sprintf('%3.3f', listNumbers(end)));   
   
    mainFolder  = rootFolder();
               
    resultsFolder = fullfile(mainFolder, 'results2Points', testType);
    if ~isdir(resultsFolder)
        mkdir(resultsFolder);
    end

    % carto file (also relative to main folder)
    path2CartoFile = fullfile('data','cartoXP','handles.mat');
    
    fileName = fullfile(resultsFolder, ['log', datestr(now), '.csv']);
    fid = fopen(fileName , 'a');

    
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


% find index for every xyz point
spatialLocation = [information.x; information.y; information.z]';

% define Euclidean distance between a point and a set of points
distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));

[~, pointAIdx] = min(distPoint2PointSet(spatialLocation, pointA));
[~, pointBIdx] = min(distPoint2PointSet(spatialLocation, pointB));
[~, pointCIdx] = min(distPoint2PointSet(spatialLocation, pointC));

allPointsIdx = [pointAIdx, pointBIdx, pointCIdx];


% Set RNG for reproducibility, keeping the old settings to be
% restored at the end. 
oldRNG = RandStream.getGlobalStream;
cleaners.RNG =...
  onCleanup(@()RandStream.setGlobalStream(oldRNG));

newRNG = RandStream('mrg32k3a', 'Seed', 0);
RandStream.setGlobalStream(newRNG);

%% where to apply fastICA.

% if isa(listOfPoints, 'char') && strcmpi(listOfPoints, 'all')
%     allPoints = unique(information.npunto);  
%     % let's assume numeric otherwise    
% else
%     % check the requested point belongs to the list. Otherwise, just wipe
%     % it out.  
%     currentPoints = unique(information.npunto);
%     presenPoints  = ismember(listOfPoints, currentPoints);
%     allPoints     = listOfPoints(presenPoints);   
%     if listOfPoints(~presenPoints)
%         fprintf('Warning: Point %d not found\n', listOfPoints(~presenPoints));
%     end
% end

numPoints      = numel(allPointsIdx);
numSelectedICs = numel(listOfSelectedIcs);
testResults    = cell(numSelectedICs, numPoints);


% Populate the header
fprintf(fid, '%s\t%s\t%s\t%s\ts%s\t%s\n', ...
    'Point', 'Input signals', 'input eKurtosis', ...
    'number of ICs', 'output eKurtosis', 'output peak Freq');




for idxPoint=1:numPoints


currentResults = [];
currentPoint   = information.npunto(allPointsIdx(idxPoint));
idx      = allPointsIdx(idxPoint);
ecgData  = data(idx);


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


%% sort data by eKurtosis to help with the representation.

kurtosisDataIn = kurtosis(selDataIn, 1, 2) - 3;
[sortedKurtDataIn, idxKurt] = sort(kurtosisDataIn) ;
selDataIn      = selDataIn(idxKurt, :);

signalNamesICA     = lisftOfSignals(idxKurt);
numInputSignalsICA = numel(signalNamesICA);

% Save some input info at this point
currentResults.point      = currentPoint;
currentResults.inputNames = signalNamesICA;
currentResults.eKurtInput = sortedKurtDataIn;

%% Preprocess here

%   TEST 1:
%
% The preprocessing consists of bandpass filtering (finite impulse response
% (FIR), 40�250 Hz, order 40, Kaiser window), rectification, and lowpass
% filtering (FIR, 0�20 Hz, order 40, Kaiser
% window)~cite{AFPreproc:richter2011novel}.
% filterOrder    = 128;
% dataSignalProc = filterFIR(selDataIn, dataIn.fs, ...
%         'filterOrder'     , filterOrder, ...
%         'filterType'      ,'high'      , ...
%         'filterFrequency' , 1.0       , ...
%         'showFilterDesign', true);

% dataForICA = filterFIR(dataSignalProc, dataIn.fs, ...
%         'filterOrder'     , filterOrder, ...
%         'filterType'      ,'low'       , ...
%         'filterFrequency' , 40         , ...        
%         'showFilterDesign', false);


% Don't use band-pass as the design is a bit agressive with the phase
% fcFilter = [0.5, 40];
% dataForICA = filterChebyShev(selDataIn, dataIn.fs, ...
%         'filterOrder' , 8, 'filterRipple', 30, ...
%         'filterFrequency' , fcFilter, ...
%         'filterType'      ,'bandpass', ...
%         'showFilterDesign', true);    
    

    
dataSignalProcs = filterChebyShev(selDataIn, dataIn.fs, ...
        'filterOrder' , 6, 'filterRipple', 30, ...
        'filterFrequency' , 3.0, ...
        'filterType'      ,'high', ...
        'showFilterDesign', false);     
 
dataForICA = filterChebyShev(dataSignalProcs, dataIn.fs, ...
        'filterOrder' , 6, 'filterRipple', 30, ...
        'filterFrequency' , 40.0 , ...
        'filterType'      ,'low', ...
        'showFilterDesign', false); 
   
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


%% run fastICA and sort output by excess Kurtosis

%covariance = cov(dataForICA');
if numFlyBys==1
    %[signalICA, A, W] = fastica(dataForICA);
[signalICA, A, W] = fastica(dataForICA, 'approach', fastICA.approach, ...
       'g', fastICA.g, 'stabilization', fastICA.stabilization); 
%  'g', fastICA.g, 'stabilization', fastICA.stabilization, ...
%  'numOfIC', 5);
   
end


if numFlyBys==2

% [~, initA, ~] = fastica(dataForICA, 'approach', fastICA.approach, ...
%     'g', fastICA.g, 'stabilization', fastICA.stabilization, 'numOfIC', 4);

[~, initA, ~] = fastica(dataForICA, 'approach', fastICA.approach, ...
    'g', fastICA.g, 'stabilization', fastICA.stabilization);

[signalICA, A, W] = fastica(dataForICA, 'initGuess', initA, ...
    'approach', fastICA.approach, ...
    'g', fastICA.g, 'stabilization', fastICA.stabilization);
    %'g', fastICA.g, 'stabilization', fastICA.stabilization, 'numOfIC', 5);
end

kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);
icaSorted = signalICA(idxKurt, :);

numOutputSignalsICA = numel(kurtosisICA);

if numInputSignalsICA ~= numOutputSignalsICA
  warning('fastICA found %d independent components (%d signals slashed)\n', ...
  numOutputSignalsICA, numInputSignalsICA- numOutputSignalsICA); %#ok<WNTAG>
end

% Save some input info at this point
currentResults.numInputSignalsICA  = numInputSignalsICA;
currentResults.kurtosisIC          = eKurtValues';
currentResults.reconsInputFreqPeak = zeros(1, numInputSignalsICA);


%% Once we have the data, iterate to get the selected ICs
for currentSelIC = 1:numSelectedICs

        icaFeatures = struct;
        selectedICs = listOfSelectedIcs{currentSelIC};

        figTitle = sprintf('filter3Hz_Point %d_ICs %s_%s', currentPoint, ...
                    commaString(selectedICs), signalsStr(lisftOfSignals));


        %% Re-construct signals with selected components
        %
        %  In this step we set to zero all of the sources regarded as VA and
        %  reconstruct observations as $hat{x} = A \times s$

        % Get the observations from the sources. 
        %   Set to zero all the ICs, and add back the ones picked.
        %   Keep the original order as they were displayed by eKurt order.
        selectedSources = zeros(size(signalICA));
        if selectedICs <= numOutputSignalsICA
            selectedSources(idxKurt(selectedICs), :) = icaSorted(selectedICs, :);
        end
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

        currentHandle = 1;

        for idx=1:numOutputSignalsICA

          outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));

          % plot time signals
          h0 = handlesArray(currentHandle); 
          plot(h0, selDataIn(idx, :), 'Color', [0 1 1], 'LineWidth', 2.5);  
          hold(h0);
          plot(h0, dataForICA(idx,:), 'LineWidth', 1.5, 'Color', [0 0 1]);
          title(h0, signalNamesICA{idx});
          axis(h0,'tight');
          plot(h0, dataForICAHat(idx,:), 'LineWidth', 1.5, 'Color', 'r');


          currentHandle = addOne(currentHandle);
          h1 = handlesArray(currentHandle);

          % Plot PSDs INPUT data
          [inputPSD, inputPSDOmega] = ...
          getPowerSpectralDensity(dataForICA(idx,:), dataIn.fs, ...
              'windowType'   , 'rectangular'            , ...
              'windowSize'   , length(dataForICA(idx,:)), ...
              'lengthFFT'    , 8192                     , ...
              'numberOverlap', []                       , ...
              'psdCutOff'    , 30);

          inputPSD = normalise(inputPSD);
          

        % hbar1 = bar(h1, inputPSDOmega, inputPSD, 'b', 'EdgeColor', [0 0 0]);    
        % children = get(hbar1, 'children');
        % alpha(children, 0.8);
        % hold(h1);
        % hbar2 = bar(h1, inputHatPSDOmega, inputHatPSD, 'r');  
        % children = get(hbar2, 'children');
        % alpha(children, 0.8);      
        % axis(h1, 'tight');


        plot(h1, inputPSDOmega, inputPSD, 'LineWidth', 2);
        [~, inputFreqPeakIdx] = max(inputPSD);

        %      hold(h1);
        %       stem(h1, inputPSDOmega(inputFreqPeakIdx), inputPSD(inputFreqPeakIdx), ...
        %           'LineWidth', 2, 'Color', [0 0.2 1]);    

        %       annotationX = inputPSDOmega(inputFreqPeakIdx)*[1,1.05];
        %       annotationY = inputPSD(inputFreqPeakIdx)*[1,-1.05];
        %       annotation(h, 'textarrow', annotationX, annotationY,'String','y = x ')

        
        % Get the PSD of the reconstructed signals
        %
        % TO-DO: Actually there's no need to do the FFT again as
        % $ FT(hat{x}) = A \times FT(s) $
        [inputHatPSD, inputHatPSDOmega] = ...
            getPowerSpectralDensity(dataForICAHat(idx,:), dataIn.fs, ...
            'windowType'   , 'rectangular', ...
            'windowSize'   , length(dataForICAHat(idx,:)), ...
            'lengthFFT'    , 8192, ...
            'numberOverlap', [], ...
            'psdCutOff'    , 30);

        inputHatPSD = normalise(inputHatPSD);
        hold(h1, 'on');
        plot(h1, inputHatPSDOmega, inputHatPSD, 'LineWidth', 2.0, 'Color', 'r');
        [~, inputHatFreqPeakIdx] = max(inputHatPSD);
        hold(h1, 'on');
        %       stem(h1, inputHatPSDOmega(inputHatFreqPeakIdx), inputHatPSD(inputHatFreqPeakIdx), ...
        %       'LineWidth', 1.5, 'Color', 'r');
        axis(h1, 'tight');

        psdText  = sprintf('[Or,Rc]=[%3.2f,%3.2f] Hz', ...
          inputPSDOmega(inputFreqPeakIdx), inputHatPSDOmega(inputHatFreqPeakIdx));
        title(h1, psdText);
        set(h1, 'xtick',[]);
        set(h1, 'xticklabel',[]);
        set(h1, 'ytick',[]);
        set(h1, 'yticklabel',[]);      
        
        currentResults.reconsInputFreqPeak(idx) = ...
            inputHatPSDOmega(inputHatFreqPeakIdx);

        currentHandle = addOne(currentHandle);
        h2 = handlesArray(currentHandle);

        currentHandle = addOne(currentHandle);
        h3 = handlesArray(currentHandle);

        [icaFeatures(idx).inputDataStats, icaFeatures(idx).powSpectDens, ...
        icaFeatures(idx).psdOmega, icaFeatures(idx).psdPeaks] = ...
        surveyIndependentComponents( icaSorted(idx, :), outSignalName, dataIn.fs, h2, h3, []);

        if any(ismember(idx, selectedICs))
          set(h2, 'LineWidth', 2, 'Color', [1 1 0.9]);
        end

        currentHandle = addOne(currentHandle);

        end

        set(gcf,'units','normalized','outerposition',[0 0 1 1]);

        % save figure
        %set(currentFig, 'units', 'normalized', 'outerposition',[0 0 1.5 1.5]);
        set(currentFig, 'Name' , figTitle);
        
        figResultsFolder = fullfile(resultsFolder, ['ics_', commaString(selectedICs)]);
        if ~isdir(figResultsFolder)
            mkdir(figResultsFolder);
        end

        figPathName = fullfile(figResultsFolder, [figTitle, '.png']);
        hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
        close(currentFig);
        fprintf('\n');


   testResults{currentSelIC, currentPoint} = currentResults;
    
   fprintf(fid, '%d\t%s\t%s\t%s\t%s\t%s\n'     , ...
   currentResults.point                  , ...
   signalsStr(currentResults.inputNames) , ...
   floatStr(currentResults.eKurtInput), ...
   commaString(selectedICs), ...
   floatStr(currentResults.kurtosisIC), ...
   floatStr(currentResults.reconsInputFreqPeak));

end
    
end
fclose(fid);
save('testResultsBoth.mat', 'testResults');