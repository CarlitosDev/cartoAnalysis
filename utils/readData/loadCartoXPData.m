function [data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile)

% LOADCARTOXPDATA from CARTO XP
% 
% Load CartoXP data using Marga's scripts. All the files are expected to
% live in FULLPATH2CARTOFILE
% 
% DATA and POINTSINFO are fetched from the .mat file
%
% MESHDATA is read using meshData = read_mesh_XP(meshFName); 
% The function will look for a 'RECONST.txt' file within the
% FULLPATH2CARTOFILE folder.
%
% CARTODATA   = read_Carto(cartoFName); 
% The function will look for a 'map.car' file within the FULLPATH2CARTOFILE folder.
%

%%
   
%     mainFolder  = rootFolder();
%     path2CartoFile = fullfile('data','cartoXP','handles.mat')
%     
% fullPath2CartoFile = fullfile(mainFolder, path2CartoFile);
% 
% %%
% 
% mainFolder  = rootFolder();
% path2CartoFile = fullfile('data','Auriculas','#2', 'Carto', 'EG.mat')
% fullPath2CartoFile = fullfile(mainFolder, path2CartoFile);


%% Define Euclidean distance between a point and a set of points
distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));
  


%% tide this up...
fileExist = @(fName) exist(fName, 'file') ~= 0;
[pathStr, ~, ~] = fileparts(fullPath2CartoFile);

% get mesh
meshFName  = fullfile(pathStr, 'RECONST.txt');
if fileExist(meshFName)
  
  % call Marga's reader
  meshData   = read_mesh_XP(meshFName);

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

% get carto file
cartoFName  = fullfile(pathStr, 'map.car');
if fileExist(cartoFName)

  cartoData   = read_Carto(cartoFName);
  [pointsId, sortIdx] = sort(cartoData.npunto);
  
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
  
else
	cartoData = [];
end

%


%% Bound mesh faces with data points
%
% Get the minimum distance between the barycenter of the facet and the
% spatial position of a data point. There are way more facets than data
% points, so they will share the same data-point.

if ~isempty(meshData) && ~isempty(cartoData)

  scaleFrom0To1 = @(x) (x-min(x))/(max(x)-min(x));
  
  closestDataPoint     = zeros(1, numFaces);
  closestPointDist     = zeros(1, numFaces);
  facetColours         = zeros(numFaces, 3);
  closestPointDistNorm = zeros(1, numFaces);
  
  
  % index colours to speed the loop up 
  coloursIndexed = zeros(max(cartoData.pointsId), 3);
  coloursIndexed(cartoData.pointsId, :) = cartoData.pointsColour;  

  for faceIDx = 1:numFaces
	[pointFacetDist, pointIdx] = ...
      min(distPoint2PointSet(pointsPosition, meshData.facesCentroid(faceIDx,:)));
    closestDataPoint(faceIDx)  = pointsId(pointIdx);
    closestPointDist(faceIDx)  = pointFacetDist;
    facetColours(faceIDx, :)   = coloursIndexed(pointsId(pointIdx), :);
  end
  
  % normalise distances from facet to data point
  [~,~,idxUniquePointFacet] = unique(closestDataPoint);
  for idx=1:numFaces
      currentIdx = idxUniquePointFacet == idx;
      closestPointDistNorm(currentIdx) = scaleFrom0To1(closestPointDist(currentIdx));
  end
  
  % update mesh data
  meshData.closestDataPoint     = closestDataPoint';
  meshData.closestPointDist     = closestPointDist';
  meshData.closestPointDistNorm = closestPointDistNorm';
  meshData.facetColours         = facetColours;
  
end


%%  

    tempVar = load(fullPath2CartoFile);
    fNames  = fieldnames(tempVar);

    % Output when reading Aitor's. I don't get what this guy's doing. Let's get
    % the data from the main struct to access all the data in the same way.
    if strcmpi(fNames, 'handles')
        information  = tempVar.handles.DB;
        tempVar.EG   = information.EG;
    end

    tempData = tempVar.ECG;
    tempPointsInfo = [tempData.ID];

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

