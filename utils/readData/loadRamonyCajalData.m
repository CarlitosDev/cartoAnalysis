function [data, pointsInfo, meshData, cartoData, cathererData] = ...
  loadRamonyCajalData(fullPath2CartoFile)
%
% LOADRAMONYCAJALDATA from CARTO data acquired at Ramon y Cajal hospital
%
% Load CartoXP data using Marga's scripts. All the files are expected to
% live in FULLPATH2CARTOFILE
%
% DATA and POINTSINFO are fetched from the .mat file. This mat file can be
% generated from scratch using 'READPENTARAYEXPORTFOLDER'.
%
% MESHDATA is read using meshData = read_mesh_C3(meshFName); The function
% will look for a '*.mesh' file within the FULLPATH2CARTOFILE folder.
%
% CARTODATA   = read_Carto(cartoFName); The function will look for a
% 'map.car' file within the FULLPATH2CARTOFILE folder.
%
% The file cartoFName contains the spatial position of the monopolar
% electrodes of the catherer branches, ie: one brach of the catherer with
% the electrodes 1-2, 3-4 will have two entries in the 'car' file, for
% electrodes 1 and 3 respectively.
%
% CATHERERDATA groups the data acquired with the catherer, the position of
% the electrodes, the type of the signal.


%%

[pathStr, ~, ~] = fileparts(fullPath2CartoFile);

%% Define Euclidean distance between a point and a set of points

distPoint2PointSet = ...
  @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));

fileExist = @(fName) exist(fName, 'file') ~= 0;

%%  Read the .mat file with the ECG signals.
% This .mat file is generated using readPentaRayExportFolder.m that needs
% to unzip the study file (nearly 16GB).

tempVar = load(fullPath2CartoFile);
fNames  = fieldnames(tempVar);

if ~strcmpi(fNames, 'ECG')
  errordlg('Can''t resolve data.EGC');
end

tempData       = tempVar.ECG;
tempPointsInfo = [tempData.ID];
numECG         = numel(tempData);

% Read the list of channels and electrodes. Make sure is compatible with
% upcoming data.
[electrodeNames, ~, nameToIdxMap] = listOfElectrodesRyC();

%% Rearrange data to match the electrode's placement
% Every data point contains 60 signals (aVF, M1, 20A_1, ...). Presumably we
% can exactly map the position of the catheter with the catherer signals.

cathererData = [];

if any(strcmpi(fieldnames(tempData), 'electrodesPos'))
  
  skipCartoFile  = true;
  
  % monopolar names  
  monopolarNames = regexp(sprintf('20A_%d\n', 1:20), '\n', 'split');
  monopolarNames(end) = [];
  
  % bipolar names
  tIdx = [1,2;2,3;3,4;5,6;6,7;7,8;9,10;10,11;11,12;13,14;...
    14,15;15,16;17,18;18,19;19,20];
  
  bipolarNames = {};
  bipolarRef   = {};
  
  for k=1:size(tIdx, 1)
    bipolarNames{k} = sprintf('20A_%d_%d', tIdx(k,1), tIdx(k,2));
    bipolarRef{k}   = sprintf('20A_%d', tIdx(k,2));
  end
  
  idxA = 0;  
  
  for idxPoint = 1:numECG
    
    ePos = tempData(idxPoint).electrodesPos;
    numElectrodes = height(ePos);
    pentaRayID    = tempData(idxPoint).ID;
    % If we've got 22, let's flush the first two as they don't record
    % signal.
    if numElectrodes == 22
      ePos(1:2, :)  = [];
      numElectrodes = 20;
    end
    
    signalNames = tempData(idxPoint).tipo_ECG;
    
    % read MONOPOLARS first
    for idxElectrode = 1:numElectrodes
      
      idxMap    = nameToIdxMap(monopolarNames{idxElectrode});
      idxSignal = find(strcmpi(signalNames, electrodeNames{idxMap}));
      
      if ~isempty(idxSignal)
        
        idxA  = idxA+1;
        
        cathererData(idxA).pentaRayElectrode = electrodeNames{idxMap, 2}; %#ok<*AGROW>
        cathererData(idxA).type              = electrodeNames{idxMap, 3};
        
        cathererData(idxA).x                 = ePos.X(idxElectrode);
        cathererData(idxA).y                 = ePos.Y(idxElectrode);
        cathererData(idxA).z                 = ePos.Z(idxElectrode);
        cathererData(idxA).time              = ePos.Time(idxElectrode);
        % 'internal' point descriptor
        cathererData(idxA).pointsId          = idxA;
        % 'external' point descriptor matching the ECG files
        cathererData(idxA).pentaRayPoint     = pentaRayID;
        cathererData(idxA).uniquePentaRayID  = idxPoint;
        
        cathererData(idxA).electrodeID       = idxElectrode;
        
        % find the current signal and copy it
        cathererData(idxA).signal            = tempData(idxPoint).data{idxSignal};
      end  
    end

    % read BIPOLARS here
    for idxElectrode = 1:numel(bipolarNames)
      
      idxMap    = nameToIdxMap(bipolarNames{idxElectrode});
      idxSignal = find(strcmpi(signalNames, electrodeNames{idxMap}));
      
      % Say the bipolar is placed at the 2nd point of the differential
      idxRemap = find(strcmpi(monopolarNames, bipolarRef{idxElectrode}));
      
      if ~isempty(idxSignal) && (idxRemap < numElectrodes)
        
        idxA  = idxA+1;
        
        cathererData(idxA).pentaRayElectrode = electrodeNames{idxMap, 2}; %#ok<*AGROW>
        cathererData(idxA).type              = electrodeNames{idxMap, 3};

        cathererData(idxA).x                 = ePos.X(idxRemap);
        cathererData(idxA).y                 = ePos.Y(idxRemap);
        cathererData(idxA).z                 = ePos.Z(idxRemap);
        cathererData(idxA).time              = ePos.Time(idxRemap);
        % 'internal' point descriptor
        cathererData(idxA).pointsId          = idxA;
        % 'external' point descriptor matching the ECG files
        cathererData(idxA).pentaRayPoint     = pentaRayID;
        cathererData(idxA).uniquePentaRayID  = idxPoint;
        
        cathererData(idxA).electrodeID       = idxElectrode;
        
        % find the current signal and copy it
        cathererData(idxA).signal            = tempData(idxPoint).data{idxSignal};
      end
    end
    
  end
  
end

%% Load the mesh file

% get mesh
meshFileList = dir(fullfile(pathStr, '*.mesh'));
meshFName    = fullfile(pathStr, meshFileList.name);

if fileExist(meshFName)
  
  % call Marga's reader
  meshData   = read_mesh_C3(meshFName);
  
  % Remove duplicated vertex
  meshVertex = [meshData.X; meshData.Y; meshData.Z]';
  meshFaces  = meshData.faces;
  numFaces   = length(meshFaces);
  [meshVertex, ~, indexMesh] = unique(meshVertex, 'rows');
  meshFaces = indexMesh(meshFaces);
  
  % apply some curvature-based smoothing
  doSmoothing = false;
  if doSmoothing
    normalisedCurvature = 1;
    numberIterations    = 1;
    lambdaSmoothing     = 1;
    neighbourInfluence  = 1;
    FV.vertices = meshVertex;
    FV.faces    = meshFaces;
    FV2 = smoothpatch(FV, normalisedCurvature, ...
      numberIterations, lambdaSmoothing, neighbourInfluence);
    
    meshFaces  = FV2.faces;
    meshVertex = FV2.vertices;
  end
  
  % update struct
  meshData.X     = meshVertex(:,1);
  meshData.Y     = meshVertex(:,2);
  meshData.Z     = meshVertex(:,3);
  meshData.faces = meshFaces;
  
  meshData.meshVertex = meshVertex;
  
  % get barycenter
  facesCentroid = zeros(numFaces, 3);
  
  for idx=1:numFaces
    vertices   = [meshVertex(meshFaces(idx, 1), :);
      meshVertex(meshFaces(idx, 2), :);
      meshVertex(meshFaces(idx, 3), :)];
    
    facesCentroid(idx, :) = mean(vertices, 1);
  end
  
  meshData.facesCentroid = facesCentroid;
  
else
  meshData = [];
end



%% Get the carto data info

% get the car carto file >> in this version is saved as .txt
cartoFileList = dir(fullfile(pathStr, '*car.txt'));
cartoFName    = fullfile(pathStr, cartoFileList.name);
if fileExist(cartoFName)
  
  cartoData = read_Carto(cartoFName);
  [pointsId, sortIdx] = sort(cartoData.npunto);
  
  % let's remove the field to avoid mistakes
  cartoData = rmfield(cartoData, 'npunto');
  
  pointsPosition = [cartoData.x; cartoData.y; cartoData.z];
  pointsPosition = pointsPosition(:, sortIdx)';
  
  cartoData.pointsId       = pointsId;
  cartoData.pointsPosition = pointsPosition;
  
  % Help plotting by getting the duple (point, colour) for every point
  numPoints = numel(cartoData.pointsId);
  newFig = figure();
  cMap   = colormap('lines');
  close(newFig);
  numColours = length(cMap);
  allColours = repmat(cMap, ceil(numPoints/numColours), 1);
  
  cartoData.pointsColour = allColours(1:numPoints, :);
  
  % re-order to match the data in cartoData
  [~, indInA, ~] = intersect(tempPointsInfo, pointsId);
  
  data = tempData(indInA);
  
  % map the varnames and the indices
  allPointsInfo = [data.ID];
  [allPointsUnique, allPointsIdx] = unique(allPointsInfo);
  pointToIdxMap = containers.Map(allPointsUnique, allPointsIdx);
  
  pointsInfo = [];
  pointsInfo.allPointsInfo = allPointsInfo;
  pointsInfo.pointToIdxMap = pointToIdxMap;
  
else
  cartoData = [];
end


%% Bound mesh faces with data points
%
% Get the minimum distance between the barycenter of the facet and the
% spatial position of a data point. There are way more facets than data
% points, so they will share the same data-point.

if ~isempty(meshData) && ~isempty(cartoData) && ~skipCartoFile
  
  [meshData, cartoData] = boundFacesWithPoints(meshData, cartoData); 
  
end



