% In this script, we select all the MONOPOLAR electrodes and then we 
% pick one of acquisition  points to visualise the meridian and parallel.

selectedElectrodeType = 'MON';

%% Load the RAW data

  dataFolder  = 'preprocessed data';
  egmFileName = 'PreProcessed_mapaFAECG.mat';
  fullPath    = fullfile(getDataFolder(), dataFolder, egmFileName);

%% No need to edit below

  a = load(fullPath);
  cathererData   = a.cathererData;
  allSignals     = {cathererData.type};
  presentSignals = unique(allSignals);

  fprintf('Available signals\n');
  for idx=1:numel(presentSignals)
    fprintf('  %s, ', presentSignals{idx});
    if mod(idx,4)==0, fprintf('\n'); end
  end
  fprintf('\n');

  idxSelSignal = strcmpi(allSignals, selectedElectrodeType);

%%

  subCartoAnalysis  = cartoPointAnalysis(idxSelSignal);
  subCathererData   = cathererData(idxSelSignal);
  currentDescriptor = botteronFrequency(idxSelSignal);

  subCartoAnalysis{1}.inputSignalName = selectedElectrodeType;

%%

  getCartoPoints = @(sData) [sData.x; sData.y; sData.z];

  for iAxis = ['x', 'y', 'z']
    cartoData.(iAxis) = [subCathererData.(iAxis)];
  end

  % Pass the data from the catherers into cartoData
  cartoData.pointsPosition = getCartoPoints(subCathererData)';
  cartoData.pointsId       = [subCathererData.pointsId];

  [meshData, cartoData]    = boundFacesWithPoints(meshData, cartoData);

%%

  plotMeshAndAllEGMAnalysisParallel(meshData, cartoData, ...
    subCartoAnalysis, currentDescriptor)