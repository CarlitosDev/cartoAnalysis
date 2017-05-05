%baseUngittedFolder = 'D:\Work\MATLABPromotionsModel\Matlab PhD';
baseUngittedFolder = '/Users/carlosAguilar/Documents/Matlab PhD/Data (ungitted)';

%baseUngittedFolder = 'D:\Work\MATLABPromotionsModel\Matlab PhD';

fullPath2SimulationFile = fullfile(baseUngittedFolder,  ...
  'SimulacionValencia', 'AF_AI.mat');

tempVar = load(fullPath2SimulationFile);

%%
getVertexIdx = @(vNumber, facetsData, facetsIdx) ...
facetsData(facetsIdx, vNumber);

getVerTexPos = @(vNumber, vertexData, facetsData, facetsIdx) ...
vertexData(getVertexIdx(vNumber, facetsData, facetsIdx), :);

getFacetPos = @(vertexData, facetsData, facetsIdx) ...
[getVerTexPos(1, vertexData, facetsData, facetsIdx); ...
getVerTexPos(2, vertexData, facetsData, facetsIdx); ...
getVerTexPos(3, vertexData, facetsData, facetsIdx)];


%%

meshData.faces = double(tempVar.faces_AI);
meshData.X = 35.0*double(tempVar.Nodos_AI(:, 1))';
meshData.Y = 35.0*double(tempVar.Nodos_AI(:, 2))';
meshData.Z = 35.0*double(tempVar.Nodos_AI(:, 3))';


%facets = meshData.faces(1:10, :);
facets = meshData.faces(1:end, :);

xVertexIdx = unique(facets(:,1));
yVertexIdx = unique(facets(:,2));
zVertexIdx = unique(facets(:,3));

x = 35.0*double(tempVar.Nodos_AI(:, 1));
y = 35.0*double(tempVar.Nodos_AI(:, 2));
z = 35.0*double(tempVar.Nodos_AI(:, 3));

currentVert = [x,y,z];
[numVertex, ~] = size(currentVert);
[numFacets, ~] = size(facets);

meshVertex = currentVert;
meshFaces = facets;

meshCentroid = mean(currentVert, 1);

%%
timeSlice = 1000;
facetColours = tempVar.egm_AI(timeSlice, :)';

%%
facetsIdx = 2;
getVertexIdx(1, facets, facetsIdx)
getVertexIdx(2, facets, facetsIdx)
getVertexIdx(3, facets, facetsIdx)

%%
figureLeft,
% meshCentroid = mean(currentVert, 1);
% trisurf(facets, x,y,z, ...
% facetColours, ...
% 'facealpha' , 0.35 )

% should be the same as

ax = newplot();

h = patch('faces',facets, ...
          'vertices',currentVert,...
          'facevertexcdata',facetColours,...
          'facecolor', 'flat', ...
          'edgecolor',get(ax,'DefaultSurfaceEdgeColor'), ...
          'parent', ax,...
          'facealpha', 0.8, ...
          'CDataMapping'   , 'scaled');
view(ax,3);
grid(ax,'on');

% notes carlos: colormap by default is set to 'jet'. As CDataMapping is set
% to scaled, it will take the mix and max values and scale jet colormap
% accordingly.

% Mathwork's help. m is the lenght of the current colormap
% index = fix((C-cmin)/(cmax-cmin)*m)+1;
% %Clamp values outside the range [1 m]
% index(index<1) = 1;
% index(index>m) = m;

%%
% figureRight,
% patch('Faces'          , facets, ...
%       'Vertices'       , currentVert , ...
%       'FaceVertexCData', facetColours, ...
%       'FaceColor'      , 'flat'      , ...
%       'facealpha'      , 0.8         , ...      
%       'CDataMapping'   , 'direct');
    

%% Heatmap colours

scaleFrom0To1    = @(x) (x-min(x))/(max(x)-min(x));
scaledDescriptor = scaleFrom0To1(facetColours);

% Set the maximum value to red and minimum to blue
heatMapArray      = zeros(numVertex, 3);
heatMapArray(:,1) = (1-scaledDescriptor);
heatMapArray(:,3) = scaledDescriptor;

% Map each facet to the closest data point
% numFaces         = length(meshData.faces);
% facetColours     = zeros(numFaces, 3);
% pointsId         = cartoData.pointsId;
% closestDataPoint = meshData.closestDataPoint;
% 
% 
% for faceIDx = 1:numFaces
%   pointIdx = closestDataPoint(faceIDx)==pointsId;
%   facetColours(faceIDx, :) = heatMapArray(pointIdx, :);
% end

figureRight,



patch('Faces'          , facets, ...
      'Vertices'       , currentVert , ...
      'FaceVertexCData', heatMapArray, ...
      'FaceColor'      , 'flat'      , ...
      'facealpha'      , 0.8         , ...      
      'CDataMapping'   , 'direct');

%%

if numVertex < 100

  xlabel('X-axis')
  ylabel('Y-axis')
  zlabel('Z-axis')

  deltaText = 0.1;

  for i=1:numFacets

  allFacetVert = getFacetPos(currentVert, facets, i);

  for idxPoint=1:3
  thisVert = allFacetVert(idxPoint, :);
  currentDirection = sign(thisVert-meshCentroid);
  textPosition = (deltaText*currentDirection)+thisVert;


  text(textPosition(1), textPosition(2), textPosition(3), ...
  sprintf('(%3.2f,%3.2f,%3.2f)', ...
  thisVert(1), thisVert(2), thisVert(3)));
  end

  end
  
end

%% get barycenter. Compare anonymous functions against for loop
