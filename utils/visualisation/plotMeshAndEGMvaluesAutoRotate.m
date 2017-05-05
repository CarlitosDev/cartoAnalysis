function plotMeshAndEGMvaluesAutoRotate( meshData, cartoData, pointDescriptor )

%% PLOTMESHANDEGMVALUESAUTOROTATE plot a CartoXP data mesh and adquisition points
% coulouring the facets according to pointDescriptor.
%
% Carlos Aguilar - April 2k16

isMatlabGT2014 = ~verLessThan('matlab','R2014b');

%% Set flags here

% Plot cartoXP recording points
plotCartoPoints = true;

% Every adquisition point is mapped to a set of facets. Enable below flag
% to plot in a different figure.
plotPointToFacetMap = false;

% Remove the mapping for those facets that are quite far from the
% adquisition point.
removeFarPoints = false;

% Overlay point number
plotPointNumber = true;

% Transparency of the mesh
faceAlpha = 0.9;

% Add interaction w the figure.
enableInteraction = true;

%%

% arrange carto data
xData = cartoData.pointsPosition(:,1);
yData = cartoData.pointsPosition(:,2);
zData = cartoData.pointsPosition(:,3);

numPoints = numel(cartoData.npunto);
assert(numel(pointDescriptor) == numPoints, ...
    'Number of descriptors doesn''t equal number of input points\n');

% arrange mesh data
currentFaces = meshData.faces;
currentVert  = [meshData.X, meshData.Y, meshData.Z];



%% Heatmap colours

scaleFrom0To1    = @(x) (x-min(x))/(max(x)-min(x));
scaledDescriptor = scaleFrom0To1(pointDescriptor);

% Set the maximum value to red and minimum to blue
heatMapArray      = zeros(numPoints, 3);
heatMapArray(:,1) = scaledDescriptor;
heatMapArray(:,3) = (1-scaledDescriptor);


% Map each facet to the closest data point
numFaces         = length(meshData.faces);
facetColours     = zeros(numFaces, 3);
pointsId         = cartoData.pointsId;
closestDataPoint = meshData.closestDataPoint;


for faceIDx = 1:numFaces
  pointIdx = closestDataPoint(faceIDx)==pointsId;
  facetColours(faceIDx, :) = heatMapArray(pointIdx, :);
end


% remove facets that lie far from data-point
if removeFarPoints
    
    meanDist = mean(meshData.closestPointDist);
    stdDist  = std(meshData.closestPointDist) ;
    idxFarPoints = meshData.closestPointDist > (meanDist+2*stdDist);
    facetColours(idxFarPoints, 1) = 1;
    facetColours(idxFarPoints, 2) = 1;
    facetColours(idxFarPoints, 3) = 1;

end

%% Plot
% Create a figure where points and geometry are overlayed. If the flag for
% interacting with the figure is being set to one, then play around with
% the figure properties and the callback. Tested in 2k14 and 2k15b.
currentFigure = figure();

% just plot carto points. No need to hold any handle at this stage.
plot3(xData, yData, zData, ...
      'o', 'MarkerSize'  , 12, ...
      'MarkerFaceColor'  , [.49 1 .63]);

hold on,
patch('Faces'          , currentFaces, ...
      'Vertices'       , currentVert , ...
      'FaceVertexCData', facetColours, ...
      'FaceColor'      , 'flat'      , ...
      'facealpha'      , faceAlpha   , ...      
      'CDataMapping'   , 'direct');
  

  
%  Let's add a customised colorbar

    [~, idxSorted] = sort(scaledDescriptor);
    maxValue       = max(pointDescriptor);
    cbarData.heatMapArray = heatMapArray(idxSorted, :); 
    cbarData.Ticks = linspace(0, maxValue, 8);
    cbarData.TickLabels = cellstr(num2str(cbarData.Ticks'));
    
    if isMatlabGT2014
        
        currentAxes = gca();
        currentAxes.UserData = cbarData;
    
        colormap(currentAxes, cbarData.heatMapArray);
        colorbar(currentAxes, 'Ticks', cbarData.Ticks./maxValue , ...
             'TickLabels', cbarData.TickLabels);
     
     end
     
%%

if plotPointToFacetMap
    figure();
    % plot carto points
    plot3(xData, yData, zData, ...
        'o', 'MarkerSize' , 12, ...
        'MarkerFaceColor',[.49 1 .63]);

    % plot mesh
    hold on,
    patch('Faces'          , currentFaces         , ...
          'Vertices'       , currentVert          , ...
          'FaceVertexCData', meshData.facetColours, ...
          'FaceColor'      , 'flat'               , ...
          'facealpha'      , 0.8                  , ...      
          'CDataMapping'   , 'direct');
end

%%

if plotPointNumber
    
    % Let's assume that the atrium is quite round. Get the centroid and
    % the direction of a vector from the centroid to every point. Then,
    % shift the text by 'delta' pointing outwards.
    deltaText      = 2.5;
    
    pointsCentroid = mean(cartoData.pointsPosition, 1);
    pointsCentered = ...
        bsxfun(@minus, cartoData.pointsPosition, pointsCentroid);
    
    vectorDirection = sign(pointsCentered);
    
    textPosition = cartoData.pointsPosition + deltaText.*vectorDirection;
       
    % TO-DO: Vectorise this...
    for i=1:numel(xData)

        hold on
        text(textPosition(i, 1), textPosition(i, 2), textPosition(i, 3), ...
            sprintf('%d', i));

    end
    
end

%% Spherical visualisation snippet.

    if ~plotPointNumber
        pointsCentroid = mean(cartoData.pointsPosition, 1);
        pointsCentered = ...
            bsxfun(@minus, cartoData.pointsPosition, pointsCentroid);
    end

    % Move along the azimuth. Go to spherical coordinates to sort the points by
    % azimuth, so we can move the camera along the horizontal axis.
    [azimuth, elevation, ~] = cart2sph(pointsCentered(:,1), ...
                                       pointsCentered(:,2), ...
                                       pointsCentered(:,3));

    [~, idxSorted] = sort(azimuth);
    [~, idxSorted] = sortrows([azimuth, elevation], [1,2]);
    
    xDataBis = xData(idxSorted, :);
    yDataBis = yData(idxSorted, :);
    zDataBis = zData(idxSorted, :);
    
    % I think the view has to revolte around the centroid
    xView = pointsCentered(idxSorted, 1);
    yView = pointsCentered(idxSorted, 2);
    zView = pointsCentered(idxSorted, 3);
      
    hold on

    for idx = 1:numPoints

        if idx>1
            plotHandle2.delete();
        end

        plotHandle2 = plot3(xDataBis(idx),yDataBis(idx),zDataBis(idx), ...
            'o', 'MarkerSize' , 30, ...
            'MarkerFaceColor' ,[0.8 0.8 0]);

       % point the camera to the point. TO-DO: Check the re-scaling
       view([xView(idx),yView(idx),zView(idx)]);

       pause(0.2)

    end



end


