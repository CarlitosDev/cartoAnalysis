% Visualise parallel and meridians. Plot the mesh an map the descriptor.
%% Set file names

% Path to CARTO file
dataSubFolder = 'RyC_Patient 2016_12_20';
dataFile      = 'ECG.mat';
cartoFile     = fullfile(getDataFolder(), dataSubFolder, dataFile);

% And analyse (BIP)OLAR or (MON)OPOLAR
currentType  = 'BIP';

%% work path out as it is ungitted

samplingFreq = 1e3;

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
idxMon    = strcmpi(cathererTable{:, 'type'}, currentType);
cathererTable = cathererTable(idxMon, :);

% % do here any aditional filtering. ie: only '20A_1_2'
% idxElectrode  = ismember(cathererTable{:, 'pentaRayElectrode'}, '20A_1_2');
% cathererTable = cathererTable(idxElectrode, :);
%   
%% Get the Botteron frequency

bottBandWidth = 2;
[botFrequencies, organisationIdx] = getBotteronAndOI(cathererTable{:, 'signal'}, ...
  samplingFreq, bottBandWidth);

% Just get the peak frequenc
fCutOff = 60;
maxFrequencies = getMaxFrequency(cathererTable{:, 'signal'}, ...
  samplingFreq, fCutOff);
plotDescriptorOverMeshV2(meshData, cathererTable, maxFrequencies);

%%  Select the point of interest. Set it manually and skip if needed. 

% Plot the atria with Botteron 
pointDescriptor  = botFrequencies;
plotDescriptorOverMeshV2(meshData, cathererTable, pointDescriptor);

% Plot the atria with OI
plotDescriptorOverMeshV2(meshData, cathererTable, organisationIdx);

%%
[selectedPointIdx, textPosition] = ...
  plotToSelectACathererPoint(meshData, cathererTable, pointDescriptor);

%% Calculate the parallels and meridians

pointsPosition = cathererTable{:, {'x', 'y', 'z'}};

[parallelPath, meridianPath, pathInfo] = ...
  getMeridianAndParallel(selectedPointIdx, pointsPosition);

%% Plot the meridian and parallels

 plotDataAnalysisParallelAndMeridianV3(cathererTable, meshData, ...
  selectedPointIdx, textPosition, parallelPath, meridianPath, pathInfo)
