% In this script, the electrode 'selectedElectrode' is selected and then we
% pick one of acquisition  points to visualise the meridian and parallel.

selectedElectrode = '20A_1_2';


%% Load the RAW data

  dataFolder  = 'preprocessed data';
  egmFileName = 'PreProcessed_mapaFAECG.mat';
  fullPath    = fullfile(getDataFolder(), dataFolder, egmFileName);

%% No need to edit below

  load(fullPath)
  cathererData   = a.cathererData;
  allSignals     = {cathererData.pentaRayElectrode};
  presentSignals = unique(allSignals);

  fprintf('Available signals\n');
  for idx=1:numel(presentSignals)
    fprintf('  %s, ', presentSignals{idx});
    if mod(idx,4)==0, fprintf('\n'); end
  end
  fprintf('\n');

  idxSelSignal = strcmpi(allSignals, selectedElectrode);

%%

  subCartoAnalysis  = cartoPointAnalysis(idxSelSignal);
  subCathererData   = cathererData(idxSelSignal);
  currentDescriptor = botteronFrequency(idxSelSignal);

  subCartoAnalysis{1}.inputSignalName = selectedElectrode;

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