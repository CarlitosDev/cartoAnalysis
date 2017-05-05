% Visualise parallel and meridians over the catheter points (not over the
% mesh).


%% Set file names

% Path to CARTO file
dataSubFolder = 'RyC_Patient 2016_12_20';
dataFile      = 'ECG.mat';
cartoFile     = fullfile(getDataFolder(), dataSubFolder, dataFile);

% And analyse (BIP)OLAR or (MON)OPOLAR
currentType  = 'BIP';

% Enable below to manually pick a point
enableInteraction = true;

% xcorr analysis
maximumLag = 400;

samplingFreq = 1e3;
%% work path out as it is ungitted

matchingFNames = which(cartoFile);
if isempty(matchingFNames)
  errordlg('Can''t resolve path to file. Run initialise()');
  return;
end
fullPath2CartoFile  = matchingFNames;
[~, srcDataName, ~] = fileparts(cartoFile);

% load loadRamonyCajalData (same interface for loadSimulation/loadCartoXPData)
[data, pointsInfo, meshData, cartoData, cathererData] = ...
  loadRamonyCajalData(fullPath2CartoFile);

cathererTable = struct2table(cathererData);

%% Select the data of interest

idxType       = strcmpi(cathererTable{:, 'type'}, currentType);
cathererTable = cathererTable(idxType, :);

% % do here any aditional filtering. ie: only '20A_1_2'
idxElectrode   = ismember(cathererTable{:, 'pentaRayElectrode'}, '20A_1_2');
cathererTable  = cathererTable(idxElectrode, :);
   
pointsPosition = cathererTable{:, {'x', 'y', 'z'}};
textPosition   = pointsPositionAsText(pointsPosition);

%% Compute a descriptor over the current points

% Get the Botteron frequency
bottBandWidth = 2;
[botFrequencies, organisationIdx] = getBotteronAndOI(cathererTable{:, 'signal'}, ...
  samplingFreq, bottBandWidth);


%% Plot the descriptor

pointDescriptor = botFrequencies;
figTitle        = 'Botteron frequencies';

% The function PLOTDESCRIPTOROVERMESHV2 plots the descriptor and maps the
% points to the current mesh.
plotDescriptorOverMeshV2(meshData, cathererTable, ...
  pointDescriptor, figTitle);

% The function PLOTDESCRIPTOROVERPOINTSINTERACTION plots the given
% descriptor over the actual acquisition points and waits until the user
% selects a point to calculate the parallel/meridians.
if enableInteraction
  selectedPointIdx = ...
    plotDescriptorOverPointsInteraction(cathererTable, ...
    pointDescriptor, figTitle);
else
  plotDescriptorOverPoints(cathererTable, pointDescriptor, figTitle); %#ok<UNRCH>
end

%% Calculate the parallels and meridians

[parallelPath, meridianPath, pathInfo] = ...
  getMeridianAndParallel(selectedPointIdx, pointsPosition);

%% Plot the meridian and parallels
doPlotMesh = false;

plotDataAnalysisParallelAndMeridianV4(cathererTable, doPlotMesh, meshData, ...
  selectedPointIdx, textPosition, parallelPath, meridianPath, pathInfo);

%% Test xcorr

% xcorr over the parallel
f    = figure();
hAxL = subplot(1,2,1);
hAxR = subplot(1,2,2);

% get the closed path
idxClosedPath = [parallelPath(1).idxPointA, parallelPath.idxPointB];

% get the signals bounded to the path
parallelBulkData = cathererTable{idxClosedPath, 'signal'};
parallelsData    = getXCorr(parallelBulkData, maximumLag);

numSamples  = size(parallelsData, 2);
timeLength  = numSamples/samplingFreq;

doScaling  = true;
pathType   = 'autocorrelation along the parallel';
pointsList = cathererTable.pointsId(pathInfo.allParallelPointsIdx);

plotSignalsIn3D(hAxL, timeLength, numSamples, ...
  parallelsData, pathType, parallelPath, pointsList, doScaling);

% xcorr over the meridian
% get the closed path
idxClosedPath = [meridianPath(1).idxPointA, meridianPath.idxPointB];

% get the signals bounded to the path
meridianBulkData = cathererTable{idxClosedPath, 'signal'};
meridiansData    = getXCorr(meridianBulkData, maximumLag);

numSamples  = size(meridiansData, 2);
timeLength  = numSamples/samplingFreq;

doScaling   = true;
pathType    = 'autocorrelation along the meridian';
pointsList  = cathererTable.pointsId(pathInfo.allMeridianPointsIdx);

plotSignalsIn3D(hAxR, timeLength, numSamples, ...
  meridiansData, pathType, meridianPath, pointsList, doScaling);

% title('Autocorrelation');

%% Test Botteron

% Parallel
[parallelsBott, parallelsBotFreqs] = ...
  getBotteronSignals(parallelBulkData, samplingFreq);

f2    = figure();
h2AxL = subplot(1,2,1);
h2AxR = subplot(1,2,2);

numSamples = size(parallelsBott, 2);
timeLength = numSamples/samplingFreq;

doScaling  = true;
pathType   = 'Botteron along the parallel';
pointsList = cathererTable.pointsId(pathInfo.allParallelPointsIdx);

plotSignalsIn3D(h2AxL, timeLength, numSamples, ...
  parallelsBott, pathType, parallelPath, pointsList, doScaling);

% Meridian
[meridiansBott, meridiansBotFreqs] = ...
  getBotteronSignals(meridianBulkData, samplingFreq);

numSamples = size(meridiansBott, 2);
timeLength = numSamples/samplingFreq;

doScaling  = true;
pathType   = 'Botteron along the meridian';
pointsList = cathererTable.pointsId(pathInfo.allMeridianPointsIdx);
  
plotSignalsIn3D(h2AxR, timeLength, numSamples, ...
  meridiansBott, pathType, meridianPath, pointsList, doScaling);

% title('Botteron signals');
