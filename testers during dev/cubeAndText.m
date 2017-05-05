%% Generate dummy data

d = [-1 1];
[x,y,z] = meshgrid(d,d,d); % a cube
x = x(:);
y = y(:);
z = z(:);

currentVert    = [x,y,z];
[numVertex, ~] = size(currentVert);

facets = ...
  [1,2,3; ...
   3,5,7; ...
   1,5,3; ...
   1,5,2; ...
  ];

[numFacets, ~] = size(facets);

numCartoPoints = 2;
idxPoints      = randi(max(numVertex), [1, numCartoPoints]);
cartoPoints    = currentVert(idxPoints, :) -1 + 2*rand(numCartoPoints, 3);

%%

figureLeft,
plot3(x,y,z, '*')

meshCentroid = mean(currentVert, 1);

hold on, 
trisurf(facets, x,y,z, ...
        'facealpha' , 0.5 )

hold on
plot3(meshCentroid(1),meshCentroid(2),meshCentroid(3), 'r*')
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')

meshCentroid = mean(currentVert, 1);

deltaText = 0.2;

for i=1:numel(x)
    
    thisVert = [x(i), y(i), z(i)];
    currentDirection = sign(thisVert-meshCentroid);
    textPosition     = (deltaText*currentDirection)+thisVert;

    text(textPosition(1), textPosition(2), textPosition(3), ...
        sprintf('Point %d', i));
    
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

end

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

% plot mesh
hold on,

patch('Faces'          , facets, ...
      'Vertices'       , currentVert , ...
      'FaceVertexCData', facetColours, ...
      'FaceColor'      , 'flat'      , ...
      'facealpha'      , 0.8         , ...      
      'CDataMapping'   , 'direct');