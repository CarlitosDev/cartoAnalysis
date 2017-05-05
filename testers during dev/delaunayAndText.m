%% Generate dummy data following trimesh example

numPointsPerDim = 5;

[x,y]  = meshgrid(1:numPointsPerDim, 1:numPointsPerDim);
z      = peaks(numPointsPerDim);
%z = z.*1e5;
facets = delaunay(x,y);

x = x(:);
y = y(:);
z = z(:);

currentVert    = [x,y,z];
[numVertex, ~] = size(currentVert);
[numFacets, ~] = size(facets);

%%
numCartoPoints = 5;
% use random points
idxPoints      = randi(max(numVertex), [1, numCartoPoints]);
cartoPoints    = currentVert(idxPoints, :) -1 + 2*rand(numCartoPoints, 3);
% use actual points
%cartoPoints = [currentVert(1, :);currentVert(24, :);currentVert(12, :);currentVert(4, :);currentVert(19, :)]

%% Plot the point number using the mesh centroid as a reference.

figureLeft,
plot3(x,y,z, '*')

meshCentroid = mean(currentVert, 1);

hold on, 
trisurf(facets, x,y,z, ...
        'facealpha' , 0.35 )

hold on
plot3(meshCentroid(1),meshCentroid(2),meshCentroid(3), 'r*')
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')

deltaText = 0.1;

for i=1:numel(x)
    
    thisVert = [x(i), y(i), z(i)];
    currentDirection = sign(thisVert-meshCentroid);
    textPosition     = (deltaText*currentDirection)+thisVert;

%     text(textPosition(1), textPosition(2), textPosition(3), ...
%         sprintf('Point %d', i));

  text(textPosition(1), textPosition(2), textPosition(3), ...
        sprintf('P-%d (%3.2f,%3.2f,%3.2f)', i, ...
        x(i), y(i), z(i)));
    
end

%% Get the barycenter
facetCentroid  = zeros(numFacets, 3); 

for idx=1:numFacets
  
  idxV1 = facets(idx, 1);
  idxV2 = facets(idx, 2);
  idxV3 = facets(idx, 3);

  v1 = currentVert(idxV1, :);
  v2 = currentVert(idxV2, :);
  v3 = currentVert(idxV3, :);

  facetPosition = [v1;v2;v3];

  facetCentroid(idx, :) = ...
    mean(facetPosition, 1);
  
%   text(facetCentroid(idx, 1), ...
%        facetCentroid(idx, 2), ...
%        facetCentroid(idx, 3), ...
%         sprintf('C-%d', idx));

cV1 = facetCentroid(idx, 1);
cV2 = facetCentroid(idx, 2);
cV3 = facetCentroid(idx, 3);

  text(cV1, cV2, cV3, ...
        sprintf('C-%d (%3.2f,%3.2f,%3.2f)', idx, ...
        cV1, cV2, cV3));


end

%% Do the same using inline methods
getVertexIdx = @(vNumber, facetsData, facetsIdx) facetsData(facetsIdx, vNumber);

getVerTexPos = @(vNumber, vertexData, facetsData, facetsIdx) ...
  vertexData(getVertexIdx(vNumber, facetsData, facetsIdx), :);

getFacetPos  = @(vertexData, facetsData, facetsIdx)    ...
  [getVerTexPos(1, vertexData, facetsData, facetsIdx); ...
   getVerTexPos(2, vertexData, facetsData, facetsIdx); ...
   getVerTexPos(3, vertexData, facetsData, facetsIdx)];

getFacetPos(currentVert, facets, idx)


%% plot the barycenter

hold on
plot3(facetCentroid(:, 1),facetCentroid(:, 2),facetCentroid(:, 3), ...
  'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g')


%% Get colourmap

% Help plotting by getting the duple (point, colour) for every point
newFig = figure(); 
cMap   = colormap('lines'); 
close(newFig);
numColours = length(cMap);
allColours = repmat(cMap, ceil(numCartoPoints/numColours), 1);

pointsColour = allColours(1:numCartoPoints, :);

%% Get the minimum distance between carto points and barycenter
distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));
  
closestDataPoint = zeros(1, numFacets);
closestPointDist = zeros(1, numFacets);
facetColours     = zeros(numFacets, 3);

% index colours to speed the loop up 
coloursIndexed = zeros(numCartoPoints, 3);
coloursIndexed(1:numCartoPoints, :) = pointsColour;  

for faceIDx = 1:numFacets
[pointFacetDist, pointIdx] = ...
    min(distPoint2PointSet(cartoPoints, facetCentroid(faceIDx,:)));
  closestDataPoint(faceIDx)  = cartoPoints(pointIdx);
  closestPointDist(faceIDx)  = pointFacetDist;
  facetColours(faceIDx, :)   = coloursIndexed(pointIdx, :);
  
end

%% Generate figure
figureRight,

% plot carto points
plot3(cartoPoints(:, 1), cartoPoints(:, 2), cartoPoints(:, 3), ...
    'o', 'MarkerSize' , 12, ...
    'MarkerFaceColor',[.49 1 .63]);

for i=1:numCartoPoints
    
    thisVert = [cartoPoints(i, 1), cartoPoints(i, 2), cartoPoints(i, 3)];
    currentDirection = sign(thisVert-meshCentroid);
    textPosition     = (deltaText*currentDirection)+thisVert;

    text(textPosition(1), textPosition(2), textPosition(3), ...
          sprintf('P-%d (%3.2f,%3.2f,%3.2f)', i, ...
          thisVert(1), thisVert(2), thisVert(3)));
    
end  
  
  
  
  
% plot mesh
hold on,

patch('Faces'          , facets, ...
      'Vertices'       , currentVert , ...
      'FaceVertexCData', facetColours, ...
      'FaceColor'      , 'flat'      , ...
      'facealpha'      , 0.8         , ...      
      'CDataMapping'   , 'direct');
    

hold on,
    
for idx=1:numFacets
  
  if mod(idx, 1) == 0
 
  text(facetCentroid(idx, 1), ...
       facetCentroid(idx, 2), ...
       facetCentroid(idx, 3), ...
        sprintf('C-%d', idx));
  end

end

%% Review some points

printPoint = @(currentPoint) fprintf('(%3.2f,%3.2f,%3.2f)', ...
  currentPoint(1), currentPoint(2), currentPoint(3));

% show current carto points
for i=1:numCartoPoints
  fprintf('\nPoint %d ->', i);
  printPoint(cartoPoints(i,:));
end


% show current facet spatial coordinates
faceIDx = 1;

idxV1 = facets(faceIDx, 1);
idxV2 = facets(faceIDx, 2);
idxV3 = facets(faceIDx, 3);

v1 = currentVert(idxV1, :);
v2 = currentVert(idxV2, :);
v3 = currentVert(idxV3, :);

facetPosition = [v1;v2;v3];
fprintf('\nFacet %d ->', faceIDx);
for i=1:3
  fprintf('\nPoint %d ->', i);
  printPoint(facetPosition(i, :));  
end

fprintf('\nFacet centroid %d ->\n', faceIDx);
printPoint(facetCentroid);

% show euclidean distances

euclidenanDist = distPoint2PointSet(cartoPoints, facetCentroid(faceIDx,:));

[pointFacetDist, pointIdx] = min(euclidenanDist);
closestDataPoint(faceIDx)  = cartoPoints(pointIdx);
closestPointDist(faceIDx)  = pointFacetDist;
facetColours(faceIDx, :)   = coloursIndexed(pointIdx, :);
  

