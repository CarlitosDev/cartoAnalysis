function plotDescriptorOverMeshV2(meshData, cathererTable, ...
  pointDescriptor, figTitle)

%% PLOTTOSELECTACATHERERPOINT plot a CartoXP data analysis
% Plot the mesh and the analysis of the adquisition points coulouring the
% facets according to pointDescriptor. Then waits for the user to pick a
% point.
%
% Carlos Aguilar - April 2k17

% v2017a has got the implicit expansion
distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));
 
%%

numPoints = height(cathererTable);

% arrange catherer data
pointsPosition = cathererTable{:, {'x', 'y', 'z'}};
xData          = cathererTable{:, 'x'};
yData          = cathererTable{:, 'y'};
zData          = cathererTable{:, 'z'};

%% Heatmap colours
maxDesc = max(pointDescriptor);
minDesc = min(pointDescriptor);

% Map each vertex to the closest data point
numVertex        = numel(meshData.X);
vertexColours    = zeros(numVertex, 1);
allVertex        = [meshData.X, meshData.Y, meshData.Z];

for vertexIDx = 1:numVertex
  [~, pointIdx] = ...
    min(distPoint2PointSet(pointsPosition, allVertex(vertexIDx,:)));
  vertexColours(vertexIDx) = pointDescriptor(pointIdx);
end

%% Plot
% Create a figure where points and geometry are overlayed. If the flag for
% interacting with the figure is being set to one, then play around with
% the figure properties and the callback. Tested in 2k14 and 2k15b.

cFigure = figureLeft();
colormap(cFigure, hot)
trimesh(meshData.faces, ...
  meshData.X, meshData.Y, meshData.Z, vertexColours);
cAx = gca();
caxis(cAx, [minDesc,maxDesc])
colorbar(cAx);
shading(gca(), 'interp');
title(figTitle);