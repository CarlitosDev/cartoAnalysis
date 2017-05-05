function plotAnalysisOverMeshInteractive( meshData, cartoData, cartoPointAnalysis, pointDescriptor )

%% PLOTANALYSISOVERMESHINTERACTIVE plot a CartoXP data mesh and adquisition points
% coulouring the facets according to pointDescriptor.
%
% Carlos Aguilar - May 2k16


%% Set flags here

% Every adquisition point is mapped to a set of facets. Enable below flag
% to plot in a different figure.
plotPointToFacetMap = false;

% Remove the mapping for those facets that are quite far from the
% adquisition point.
removeFarPoints = false;

% Overlay point number
plotPointNumber = true;

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
currentVert  = [meshData.X,meshData.Y,meshData.Z];

% save the id of the points to add them to the plot
allPointsId  = zeros(1, numPoints);

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

if enableInteraction
  
  currentAxes = gca;
  plotHandle  = [];
  
  % plot carto points one by one
  for pointIdx = 1:numel(xData)
    
    % plot independently so we keep the handle
    plotHandle(pointIdx) = ...
      plot3(currentAxes, xData(pointIdx), yData(pointIdx), zData(pointIdx), ...
      'o', 'MarkerSize' , 12, ...
      'MarkerFaceColor' ,[.49 1 .63]);
  
    hold on,
    
    setappdata(plotHandle(pointIdx), 'cartoAnalysis', cartoPointAnalysis{pointIdx});    
%     
%     set(plotHandle(pointIdx), ...
%         'HitTest', 'on', ...
%         'ButtonDownFcn', @(~,~) plotCartoAnalysisCallback(plotHandle(pointIdx)));
    set(plotHandle(pointIdx), ...
        'HitTest', 'on', ...
        'ButtonDownFcn', @(~,~) plotCartoAnalysisOverlayCallback(plotHandle(pointIdx)));      
    
  end
   
else

    % just plot carto points. No need to hold any handle at this stage.
    plot3(xData, yData, zData, ...
          'o', 'MarkerSize'  , 12, ...
          'MarkerFaceColor'  , [.49 1 .63]);

end
  


% plot mesh
hold on,
patch('Faces'          , currentFaces, ...
      'Vertices'       , currentVert , ...
      'FaceVertexCData', facetColours, ...
      'FaceColor'      , 'flat'      , ...
      'facealpha'      , 1.0         , ...      
      'CDataMapping'   , 'direct');
  

  
%  Use this code to add a colorbar
    [~, idxSorted] = sort(scaledDescriptor);
    maxValue       = max(pointDescriptor);
    cbarData.heatMapArray = heatMapArray(idxSorted, :); 
    cbarData.Ticks = linspace(0, maxValue, 8);
    cbarData.TickLabels = cellstr(num2str(cbarData.Ticks'));
    
    if ~verLessThan('matlab','R2014b')
        
        currentAxes = gca();
        currentAxes.UserData = cbarData;
    
        colormap(currentAxes, cbarData.heatMapArray);
        colorbar(currentAxes, 'Ticks', cbarData.Ticks./maxValue , ...
             'TickLabels', cbarData.TickLabels);
     
     end




if plotPointNumber
    
    % Let's assume that the atrium is quite round and get the centroid and
    % get the direction of a vector from the centroid to every point. Then,
    % shift the text by 'delta' pointing outwards.
    deltaText      = 3.5;
    pointsCentroid = mean(cartoData.pointsPosition, 1);
    
    vectorDirection = ...
        sign(bsxfun(@minus, cartoData.pointsPosition, pointsCentroid));
    
    textPosition = cartoData.pointsPosition + deltaText.*vectorDirection;
       
    
    % TO-DO: Vectorise this...
    for i=1:numel(xData)
        
        allPointsId(i) = cartoPointAnalysis{i}.pointId;

        hold on
        text(textPosition(i, 1), textPosition(i, 2), textPosition(i, 3), ...
            sprintf('%d', allPointsId(i)));

    end
    
    
end

alignFigure(currentFigure, 'left');

%% Points to facets.
% Plot the mapping from every adquisition point to the facets in a new
% figure.

if plotPointToFacetMap
    figure()
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




end

