baseUngittedFolder = '/Users/carlosAguilar/Documents/Matlab PhD/Data (ungitted)';

%baseUngittedFolder = 'D:\Work\MATLABPromotionsModel\Matlab PhD';

fullPath2SimulationFile = fullfile(baseUngittedFolder,  ...
  'SimulacionValencia', 'AF_AI.mat');



tempVar = load(fullPath2SimulationFile);


  meshData.faces = double(tempVar.faces_AI);
  meshData.X     = 35.0*double(tempVar.Nodos_AI(:, 1))';
  meshData.Y     = 35.0*double(tempVar.Nodos_AI(:, 2))';
  meshData.Z     = 35.0*double(tempVar.Nodos_AI(:, 3))';


facets = meshData.faces(1:10, :);

xVertexIdx = unique(facets(:,1));
yVertexIdx = unique(facets(:,2));
zVertexIdx = unique(facets(:,3));

x    = 35.0*double(tempVar.Nodos_AI(:, 1));
y   = 35.0*double(tempVar.Nodos_AI(:, 2));
z    = 35.0*double(tempVar.Nodos_AI(:, 3));

currentVert    = [x,y,z];
[numVertex, ~] = size(currentVert);
[numFacets, ~] = size(facets);

%%
getVertexIdx = @(vNumber, facetsData, facetsIdx) ...
                 facetsData(facetsIdx, vNumber);

getVerTexPos = @(vNumber, vertexData, facetsData, facetsIdx) ...
                vertexData(getVertexIdx(vNumber, facetsData, facetsIdx), :);
              
getFacetPos  = @(vertexData, facetsData, facetsIdx)    ...
                [getVerTexPos(1, vertexData, facetsData, facetsIdx); ...
                 getVerTexPos(2, vertexData, facetsData, facetsIdx); ...
                 getVerTexPos(3, vertexData, facetsData, facetsIdx)];

%%

figureLeft,
meshCentroid = mean(currentVert, 1);
trisurf(facets, x,y,z, ...
        'facealpha' , 0.35 )
      
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')

deltaText = 0.1;

for i=1:numFacets
    
    allFacetVert = getFacetPos(currentVert, facets, i);
    
    for idxPoint=1:3
    thisVert         = allFacetVert(idxPoint, :);
    currentDirection = sign(thisVert-meshCentroid);
    textPosition     = (deltaText*currentDirection)+thisVert;


    text(textPosition(1), textPosition(2), textPosition(3), ...
          sprintf('(%3.2f,%3.2f,%3.2f)', ...
          thisVert(1), thisVert(2), thisVert(3)));
    end
    
end