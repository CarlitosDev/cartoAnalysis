function [data, pointsInfo, meshData, cartoData] = loadSimulationData(fullPath2SimFile)

% LOADSIMULATIONDATA load a mat file containing an electrical simulation of
% the atria. Load a file that contains the simulated egm data and the mesh
% where the data is recorded.
%
%
% Carlos Aguilar - February 2k16

%%
% 
% Load simulation
% baseUngittedFolder = ...
%   '/Users/carlosAguilar/Documents/Matlab PhD/Data (ungitted)';
% 
% fullPath2SimFile = fullfile(baseUngittedFolder,'SimulacionValencia', 'AF_AI.mat');
%    fullPath2SimFile = ...
%        'D:\Work\MATLABPromotionsModel\Matlab PhD\Simulacion Valencia\AF_AI.mat';
% 
% %%
% 
% mainFolder  = rootFolder();
% path2CartoFile = fullfile('data','Auriculas','#2', 'Carto', 'EG.mat')
% fullPath2CartoFile = fullfile(mainFolder, path2CartoFile);


decimationFactor = [];
% between 0 and 1.0
%decimationFactor = 0.0065;
decimationFactor = 0.005;


%% Define Euclidean distance between a point and a set of points
distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));
  
%% Do the same to get the position of a facet
getVertexIdx = @(vNumber, facetsData, facetsIdx) ...
                 facetsData(facetsIdx, vNumber);

getVerTexPos = @(vNumber, vertexData, facetsData, facetsIdx) ...
                vertexData(getVertexIdx(vNumber, facetsData, facetsIdx), :);

getFacetPos  = @(vertexData, facetsData, facetsIdx)    ...
                [getVerTexPos(1, vertexData, facetsData, facetsIdx); ...
                 getVerTexPos(2, vertexData, facetsData, facetsIdx); ...
                 getVerTexPos(3, vertexData, facetsData, facetsIdx)];

% getFacetPos(currentVert, facets, idx)
  
%%  

tempVar = load(fullPath2SimFile);
fNames  = fieldnames(tempVar);

%% Read mesh

% Current dimensions in mm
%   X => (-4.24,4.12)[8.36]  	
%   Y => (-0.08,4.81)[4.89]  
%   Z => (-18.84,-13.91)[4.93]
% 
% Normal left atria dimensions:
%   Diameter: 28-40 mm
%   Major axis: 41-61
% web.stanford.edu/group/ccm_echocardio/cgi-bin/mediawiki/index.php/Left_atrium_dimensions
%   let's scale by 20.0

  meshData.faces = double(tempVar.faces_AI);
  meshData.X     = 35.0*double(tempVar.Nodos_AI(:, 1))';
  meshData.Y     = 35.0*double(tempVar.Nodos_AI(:, 2))';
  meshData.Z     = 35.0*double(tempVar.Nodos_AI(:, 3))';

  % Arrange data
  meshVertex = [meshData.X; meshData.Y; meshData.Z]';
  meshFaces  = meshData.faces;
  numFaces   = length(meshFaces);
  
  % Remove duplicated vertex
%   [meshVertex, ~, indexMesh] = unique(meshVertex, 'rows');
%   meshFaces = indexMesh(meshFaces);
  
  % apply some curvature-based smoothing
  doSmoothing = false;
  if doSmoothing
    normalisedCurvature = 1; %#ok<UNRCH>
    numberIterations    = 6;
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
  meshData.faces      = meshFaces;
  meshData.X          = meshVertex(:,1);
  meshData.Y          = meshVertex(:,2);
  meshData.Z          = meshVertex(:,3);
  meshData.meshVertex = meshVertex;
  
  % get barycenter
  % TO-DO: No need to use a for loop in here...
  facesCentroid = zeros(numFaces, 3);

  for idx=1:numFaces
    vertices   = [meshVertex(meshFaces(idx, 1), :);
                  meshVertex(meshFaces(idx, 2), :);
                  meshVertex(meshFaces(idx, 3), :)];              

    facesCentroid(idx, :) = mean(vertices, 1);
  end

% %That's cool but way slower than the for-loop
%   for idx=1:numFaces
%     vertices = getFacetPos(meshVertex, meshFaces, idx);
%     facesCentroid(idx, :) = mean(vertices, 1);
%   end


  meshData.facesCentroid = facesCentroid;

%% Load cartoXP data (simulated EGM)

% Tricky bit - As the data is coming from a simulation, we have an EGM
% signal per vertex, so let's return all of them.    

  tempData           = tempVar.egm_AI;
  [~, numPoints] = size(tempData);

  if isempty(decimationFactor)
    % just copy over the data from the simulation   
    cartoData.pointsId       = 1:numPoints;
    cartoData.pointsPosition = meshData.meshVertex;
    cartoData.npunto         = 1:numPoints;
  else
    
    % save the source points  
    cartoData.pointsId_src       = 1:numPoints;
    cartoData.pointsPosition_src = meshData.meshVertex;
    cartoData.npunto_src         = 1:numPoints;
    cartoData.decimationFactor   = decimationFactor;
    cartoData.egmData            = tempData;
    
    % decimate the data points to fake an actual EP lab session
    numDecimatedPoints = floor(decimationFactor*numPoints);
    cartoPointsIdx     = sort(randperm(numPoints, numDecimatedPoints));
    cartoData.pointsId = cartoPointsIdx;
    cartoData.pointsPosition = meshData.meshVertex(cartoPointsIdx, :);     
    % overwrite data variables
    numPoints = numDecimatedPoints;
    tempData = tempData(:, cartoPointsIdx);
    
    cartoData.npunto = 1:numPoints;
  end
  
    pointsInfo = cartoData.pointsId;
  

  % Set the data in the same format as the actual cartoXP data
  %
  % tempData.ID
  % tempData.tipo_ecg %cell 1xN string
  % tempData.signal % cell 1xN with (N) arrays of 1xSAMPLES
  
  data = [];
  % me and cell arrays...
  for idxPoint=1:numPoints
    data(idxPoint).tipo_ECG = 'egm_sim';
    data(idxPoint).ID = pointsInfo((idxPoint));
    data(idxPoint).signal = {tempData(:, idxPoint)'};
  end
  
  
  % Help plotting by getting the duple (point, colour) for every point;
  newFig = figure(); 
  cMap   = colormap('lines'); 
  close(newFig);
  numColours = length(cMap);
  allColours = repmat(cMap, ceil(numPoints/numColours), 1);

  cartoData.pointsColour = allColours(1:numPoints, :);

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
      min(distPoint2PointSet(cartoData.pointsPosition, meshData.facesCentroid(faceIDx,:)));
    closestDataPoint(faceIDx)  = cartoData.pointsId(pointIdx);
    closestPointDist(faceIDx)  = pointFacetDist;
    facetColours(faceIDx, :)   = coloursIndexed(cartoData.pointsId(pointIdx), :);
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
