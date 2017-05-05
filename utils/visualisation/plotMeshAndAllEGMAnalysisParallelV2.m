function plotMeshAndAllEGMAnalysisParallelV2( meshData, cartoData, cartoPointAnalysis, pointDescriptor)

%% PLOTMESHANDALLEGMANALYSIS plot a CartoXP data analysis
% Plot the mesh and the analysis of the adquisition points coulouring the
% facets according to pointDescriptor.
%
% Carlos Aguilar - June 2k16

isMatlabGT2014 = ~verLessThan('matlab','R2014b');
distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));
  
%% Set flags here
% Plot the mesh as a wireframe
plotWireframe = false;

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

% Add interaction w/ the figure.
enableInteraction = true;

%%

% arrange carto data
xData = cartoData.pointsPosition(:,1);
yData = cartoData.pointsPosition(:,2);
zData = cartoData.pointsPosition(:,3);


pointsId  = cartoData.pointsId;
numPoints = numel(pointsId);
assert(numel(pointDescriptor) == numPoints, ...
    'Number of descriptors doesn''t equal number of input points\n');

% arrange mesh data
currentFaces = meshData.faces;
currentVert  = [meshData.X, meshData.Y, meshData.Z];


%% plot wireframe with EGM points

if plotWireframe
    figure();
    % plot carto points
    plot3(xData, yData, zData, ...
        'o', 'MarkerSize' , 12, ...
        'MarkerFaceColor',[.49 1 .63]);
    % plot mesh
    hold on,
    trimesh(currentFaces, ...
            meshData.X, meshData.Y, meshData.Z, ...
            'FaceColor', 'none', ...
            'EdgeColor', [0 0 1]);
end

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

currentFigure = figureLeft();

if enableInteraction
  
  currentAxes = gca;
  plotHandle  = [];
  
  % plot carto points one by one
  for pointIdx = 1:numel(xData)
    
    % plot independently so we keep the handle
    plotHandle(pointIdx) = ...
      plot3(currentAxes, xData(pointIdx), yData(pointIdx), zData(pointIdx), ...
      'o', 'MarkerSize' , 12, ...
      'MarkerFaceColor' ,[.49 1 .63]); %#ok<*AGROW>
  
    hold on,
    
    set(plotHandle(pointIdx), ...
        'HitTest', 'on', ...
        'ButtonDownFcn', @(~,~) addPointIndex(currentFigure, pointIdx));    
  end
   
else

    % just plot carto points. No need to hold any handle at this stage.
    plot3(xData, yData, zData, ...
          'o', 'MarkerSize'  , 12, ...
          'MarkerFaceColor'  , [.49 1 .63]);

end

% can't do shading as I am using direct colour mapping

patch('Faces'          , currentFaces, ...
      'Vertices'       , currentVert , ...
      'FaceVertexCData', facetColours, ...
      'FaceColor'      , 'flat'      , ...
      'facealpha'      , faceAlpha   , ...      
      'CDataMapping'   , 'direct');

  
% lighting gouraud
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
            sprintf('%d', pointsId(i)));

    end
    
end

%% Spherical visualisation snippet.

    % pause Matlab until 'continue' is pressed
    h = uicontrol('String', 'Continue', 'Position', [20 20 100 30], ...
    'Callback', 'set(gcbf, ''Name'', sprintf(''Point %d '', getappdata(gcbf,''pointIndices'')))');
    waitfor(currentFigure, 'Name');
    selectedPointIdx = getappdata(currentFigure, 'pointIndices');

    if ~plotPointNumber
        pointsCentroid = mean(cartoData.pointsPosition, 1);
        pointsCentered = ...
            bsxfun(@minus, cartoData.pointsPosition, pointsCentroid);
    end

    % Move along the azimuth. Go to spherical coordinates to sort the points by
    % azimuth, so we can move the camera along the horizontal axis.
    [currentAzimuth, currentElevation, currentRadius] = ...
      cart2sph(pointsCentered(:,1), ...
      pointsCentered(:,2), ...
      pointsCentered(:,3));
                                                   
    currentAzimuth   = currentAzimuth   + pi;
    currentElevation = currentElevation + 0.5*pi;

    % let's say we pick a point 
    
%save('fullyCrap.mat')    
    
    selectedPointAz = currentAzimuth(selectedPointIdx);
    selectedPointEl = currentElevation(selectedPointIdx);

    azimuthRange   = 10*pi/180;
    elevationRange =  5*pi/180;

    maxAzimuth     = selectedPointAz + azimuthRange;
    minAzimuth     = selectedPointAz - azimuthRange;

    maxElevation   = selectedPointEl + elevationRange;
    minElevation   = selectedPointEl - elevationRange;

    % pick every elevation point within the azimuth (include the other side!)
    validAzimuth   = currentAzimuth <  maxAzimuth & ...
                     currentAzimuth >  minAzimuth;

    % dark side of the moon...
    if maxAzimuth >= pi
        darkSideMaxAzimuth = maxAzimuth - pi;
        darkSideMinAzimuth = minAzimuth - pi;
    else
        darkSideMaxAzimuth = maxAzimuth + pi;
        darkSideMinAzimuth = minAzimuth + pi;
    end
    
    validAzimuth   = validAzimuth   | ... 
                     currentAzimuth < darkSideMaxAzimuth & ...
                     currentAzimuth > darkSideMinAzimuth;             


    validElevation = currentElevation < maxElevation & ...
                     currentElevation > minElevation;

    idxClosePoints = find(validAzimuth | validElevation);
  
    [~, idxSortedByAzEl] = ...
      sortrows([currentAzimuth(idxClosePoints), currentElevation(idxClosePoints)], [1,2]);

    idxSorted = idxClosePoints(idxSortedByAzEl);
    
    % find previous index in the re-arranged list
    sortSelectedPointIdx = find(idxSorted == selectedPointIdx);
  
    maxRadius     = max(currentRadius);
    fullElevation = 0.5*(-pi:2*pi/100:pi);
    fullAzimuth   = -pi:2*pi/100:pi;

%% create meridian rings (centered)

    [xTemp, yTemp, zTemp] = sph2cart(minAzimuth-pi, fullElevation, maxRadius);
    meridianLeftA  = [xTemp; yTemp; zTemp];
    [xTemp, yTemp, zTemp] = sph2cart(minAzimuth   , fullElevation, maxRadius);
    meridianLeftB  = [xTemp; yTemp; zTemp];
    meridianLeft   = [meridianLeftA, meridianLeftB];
    meridianLeft   = ...
        bsxfun(@plus, meridianLeft, pointsCentroid');

    % dark side of the moon...
    % if maxAzimuth >= pi
    % end

    [xTemp, yTemp, zTemp] = sph2cart(maxAzimuth-pi, fullElevation, maxRadius);
    meridianRightA = [xTemp; yTemp; zTemp];
    [xTemp, yTemp, zTemp] = sph2cart(maxAzimuth   , fullElevation, maxRadius);
    meridianRightB = [xTemp; yTemp; zTemp];
    meridianRight  = [meridianRightA, meridianRightB];
    meridianRight   = ...
        bsxfun(@plus, meridianRight, pointsCentroid');
    
%% create parallel rings

% x = r .* cos(elevation) .* cos(azimuth)
% y = r .* cos(elevation) .* sin(azimuth)
% z = r .* sin(elevation)

% that will make 'z' as a single number...

[xTemp, yTemp, zTemp] = sph2cart(fullAzimuth, minElevation-0.5*pi, maxRadius);
zTemp = zTemp.*ones(size(yTemp));
parallelDown = [xTemp; yTemp; zTemp];

parallelDown   = ...
    bsxfun(@plus, parallelDown, pointsCentroid');


[xTemp, yTemp, zTemp] = sph2cart(fullAzimuth, maxElevation-0.5*pi, maxRadius);
zTemp = zTemp.*ones(size(yTemp));
parallelUp = [xTemp; yTemp; zTemp];

parallelUp   = ...
    bsxfun(@plus, parallelUp, pointsCentroid');


%%

hold on,
plot3(meridianLeft(1,:), meridianLeft(2,:), meridianLeft(3,:));

hold on,
plot3(meridianRight(1,:), meridianRight(2,:), meridianRight(3,:));

hold on,
plot3(parallelDown(1,:), parallelDown(2,:), parallelDown(3,:));

hold on,
plot3(parallelUp(1,:), parallelUp(2,:), parallelUp(3,:));
    
  %%  Plot a new figure
    
  xDataBis = xData(idxSorted, :);
  yDataBis = yData(idxSorted, :);
  zDataBis = zData(idxSorted, :);
              
  selectedPointsData              = [];
  selectedPointsData.selectedIdx  = sortSelectedPointIdx;
  selectedPointsData.xData        = xDataBis;
  selectedPointsData.yData        = yDataBis;
  selectedPointsData.zData        = zDataBis;
  selectedPointsData.textPosition = textPosition(idxSorted, :);
  selectedPointsData.pointsId     = pointsId(idxSorted);

%%    
    
plotCartoDataAnalysisSliced(idxSorted, cartoPointAnalysis, meshData, selectedPointsData);




end




