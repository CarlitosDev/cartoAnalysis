function [meshData, cartoData] = boundFacesWithPoints(meshData, cartoData)
% Bound mesh faces with data points
%
% Get the minimum distance between the barycenter of the facet and the
% spatial position of a data point. There are way more facets than data
% points, so they will share the same data-point.

if ~isempty(meshData) && ~isempty(cartoData)
  
  scaleFrom0To1 = @(x) (x-min(x))/(max(x)-min(x));
  distPoint2PointSet = ...
  @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));

  pointsPosition = cartoData.pointsPosition;
  
  % Help plotting by getting the duple (point, colour) for every point
  numPoints = numel(cartoData.pointsId);
  newFig = figure();
  cMap   = colormap('lines');
  close(newFig);
  numColours = length(cMap);
  allColours = repmat(cMap, ceil(numPoints/numColours), 1);
  
  cartoData.pointsColour = allColours(1:numPoints, :);
  
  numFaces             = length(meshData.faces);
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
    closestDataPoint(faceIDx)  = cartoData.pointsId(pointIdx);
    closestPointDist(faceIDx)  = pointFacetDist;
    facetColours(faceIDx, :)   = coloursIndexed(cartoData.pointsId(pointIdx), :);
  end
  
  % normalise distances from facet to data point
  [~,~,idxUniquePointFacet] = unique(closestDataPoint);
  for idx=1:numFaces
    currentIdx = idxUniquePointFacet == idx;
    closestPointDistNorm(currentIdx) = ...
      scaleFrom0To1(closestPointDist(currentIdx));
  end
  
  % update mesh data
  meshData.closestDataPoint     = closestDataPoint';
  meshData.closestPointDist     = closestPointDist';
  meshData.closestPointDistNorm = closestPointDistNorm';
  meshData.facetColours         = facetColours;
  
end
