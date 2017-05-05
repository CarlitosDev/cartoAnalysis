% Discuss with JL why the data is so different when we focus on a single
% point and we visualise of the signals acquired w the catherer for that
% point.

%% Set file names

% Path to CARTO file
dataSubFolder = 'RyC_Patient 2016_12_20';
dataFile      = 'ECG.mat';
cartoFile     = fullfile(getDataFolder(), dataSubFolder, dataFile);


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

cathTable = struct2table(cathererData);

%% Play around in here

% Get a list of the points recorded with the catherer
listOfPoints = unique(cathTable{:, 'pentaRayPoint'});
% select a random one
currentPoint = listOfPoints(9);
% And analyse (BIP)OLAR or (MON)OPOLAR
currentType  = 'MON';

plotCatheterSignalsFrequencyV2(cathTable, currentPoint, currentType);

