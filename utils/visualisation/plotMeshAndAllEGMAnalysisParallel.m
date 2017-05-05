function plotMeshAndAllEGMAnalysisParallel(meshData, cartoData, cartoPointAnalysis, pointDescriptor)

%% PLOTMESHANDALLEGMANALYSIS plot a CartoXP data analysis
% Plot the mesh and the analysis of the adquisition points coulouring the
% facets according to pointDescriptor. Then waits for the user to pick a
% point where to generate the meridian/parallel analysis.
%
% Carlos Aguilar - June 2k16

isMatlabGT2014 = ~verLessThan('matlab','R2014b');
isMatlabGT2015 = ~verLessThan('matlab','R2015b');
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
faceAlpha = 0.2;

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
        'ButtonDownFcn', @(~,~) addReferencePoint(currentFigure, pointIdx));  
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
    if isMatlabGT2015
      h = uicontrol('String', 'Continue', 'Position', [20 20 100 30], ...
      'Callback', 'set(gcbf, ''Name'', sprintf(''Point %d '', getappdata(gcbf,''pointIndices'')))');
      waitfor(currentFigure, 'Name');
      selectedPointIdx = getappdata(currentFigure, 'pointIndices');
    else
      currentVer = ver('matlab');
      warningMessage = sprintf('%s %s %s %s %s\n%s\n', currentVer.Name, ...
        currentVer.Version, currentVer.Release, ...
      'does not support the current functionality.', ...
      'Input point number using the prompt.');
      title(gca, warningMessage);
      warning(warningMessage); %#ok<SPWRN>
      selectedPointIdx = input('Type a point to run the visualisation\n');
    end
    
    
    currentFigure.delete();
    
    % center the dataset
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

    
  %  Plot a new figure
    [parallelPath, meridianPath, pathInfo] = ...
        getMeridianAndParallelV4(currentAzimuth, currentElevation, ...
        currentRadius, selectedPointIdx, pointsCentered);

    % Arrange selected data
    selectedPointsData = [];
    
samplingFrequency = 1e3;
lengthSignal      = 2500;
selectedPoint     = pointsId(selectedPointIdx);
    
    % add some info
      selectedPointsData.testType        = 'manual input';
      selectedPointsData.inputSignalName = cartoPointAnalysis{1}.inputSignalName;
      selectedPointsData.fs              = samplingFrequency;
      selectedPointsData.lengthSignal    = lengthSignal;
      
    % push the data from the analysis
      selectedPointsData.selectedPoint    = selectedPoint;
      selectedPointsData.selectedPointIdx = selectedPointIdx;
      selectedPointsData.selectedAz       = currentAzimuth(selectedPointIdx);
      selectedPointsData.selectedEl       = currentElevation(selectedPointIdx);
      selectedPointsData.textPosition     = textPosition;
      selectedPointsData.pointsId         = pointsId;      
      %selectedPointsData.pointsId         = cartoData.pointsId;

      
    % add summaries of the paths
      selectedPointsData.numParallelPoints    = pathInfo.numParallelPoints;
      selectedPointsData.allParallelPointsIdx = pathInfo.allParallelPointsIdx;      
      selectedPointsData.numMeridianPoints    = pathInfo.numMeridianPoints;
      selectedPointsData.allMeridianPointsIdx = pathInfo.allMeridianPointsIdx;

    % add the paths
      selectedPointsData.parallelPath         = parallelPath;
      selectedPointsData.meridianPath         = meridianPath;
      selectedPointsData.pathInfo             = pathInfo;




%% call the plotter
    
plotDataAnalysisParallelAndMeridianV2(cartoPointAnalysis, meshData, ...
  xData, yData , zData, selectedPointsData);

end




