function  [parallelPath, meridianPath] = getMeridianAndParallelMercator(currentAzimuth, currentElevation, selectedPointIdx)
% GETMERIDIANANDPARALLELMERCATOR find paths along the azimuth and elevation
% of a (spherical) dataset for a given point making use of Mercator
% projections.
% 
% To make it easier, the spherical coordinates are projected using normal
% and transverse Mercator (azimuth and elevation respectively) to a square
% canvas so the path can be are searched in this projection. Therefore the
% distance is just an indicator of the 'closeness' between points, not an
% actual distance.
%
%
% Carlos Aguilar - 4th September 2k16

%% set up

  validDistance   = @(d) ~(isnan(d) || isinf(d));
  beVerbose       = false;
  plotPaths       = false;%true;%
  plotProjections = false;  
  debugMode       = false;
  
  numPoints  = numel(currentAzimuth);
  
  % center the projections on the point of interest
  lambdaZero = currentAzimuth(selectedPointIdx);
  phiZero    = currentElevation(selectedPointIdx);
  
%% normal mercator projection

xNormMercator = currentAzimuth - lambdaZero;
yNormMercator = log(tan(pi/4 + currentElevation/2));


%% azimuth >> Normal Mercator

  selectedPointAz  = xNormMercator(selectedPointIdx);
  selectedPointEl  = yNormMercator(selectedPointIdx);

  % Recenter from 0 to 2pi 
  idxShift = xNormMercator < selectedPointAz;
  xNormMercator(idxShift) = xNormMercator(idxShift)+(2*pi);

  % prioritise x-axis in the distance calculation
  yNormMercator = 1.5.* yNormMercator;
  
  % push the selected point at the end to close the loop
  xNormMercator(numPoints+1) = selectedPointAz + 2*pi;
  yNormMercator(numPoints+1) = selectedPointEl;


if plotProjections
    figure;
    scatter(xNormMercator, yNormMercator);
    hold on,
    scatter(xNormMercator(selectedPointIdx), ...
            yNormMercator(selectedPointIdx), ...
            'r', 'filled', 'LineWidth', 4 );
    title('Normal Mercator projection');
end  
  
%% Find the points along the azimuth
  
  parallelPath    = [];
  pointNumber     = 0;

  remainingPoints = horzcat(xNormMercator, yNormMercator);
  currentPointIdx = selectedPointIdx;
  cDist = 0;
  
  while validDistance(cDist)
  
    currentPoint    = remainingPoints(currentPointIdx, :);
    remainingPoints(currentPointIdx, :) = inf;
    
    remainingIdx    = remainingPoints(:, 1) > currentPoint(1);
    remainingPoints(~remainingIdx, :) = inf;
    
    pointNumber   = pointNumber+1;
    [cDist, cIdx] = pdist2(remainingPoints, currentPoint, 'euclidean', 'Smallest', 1);
    
    parallelPath(pointNumber).idxPointA = currentPointIdx;
    parallelPath(pointNumber).idxPointB = cIdx;
    parallelPath(pointNumber).distAB    = cDist;
    
    if beVerbose
      fprintf('Point %d->%d dist %3.2f...\n', currentPointIdx, cIdx, cDist);
    end

    % plot the walk
    if debugMode
      hold on,
      scatter(xNormMercator(cIdx), yNormMercator(cIdx), ...
              'g', 'filled', 'LineWidth', 6 );
      pause();
    end
    
    % update
    currentPointIdx = cIdx;
    
  end
    
  % If all worked out as planned, last record should be pulled out
  if parallelPath(end).idxPointA == numPoints+1
    parallelPath(end) = [];
  end

%% transverse mercator projection 

  xTransMercator = atanh(cos(currentElevation).*sin(currentAzimuth - lambdaZero));
  yTransMercator = atan(tan(currentElevation)./cos(currentAzimuth - lambdaZero)) - phiZero;
  
%% elevation >> transverse mercator projection
  
  selectedPointAz  = xTransMercator(selectedPointIdx);
  selectedPointEl  = yTransMercator(selectedPointIdx);
  
  % Recenter the poles...that's a bit of a hack
  idxShift     = yTransMercator < selectedPointEl;
  maxElevation = max(yTransMercator);
  yTransMercator(idxShift) = maxElevation - yTransMercator(idxShift);  

  % prioritise x-axis in the distance calculation
  xTransMercator = 1.5.* xTransMercator;
  
  % push the selected point at the end to close the loop
  xTransMercator(numPoints+1) = selectedPointAz;
  yTransMercator(numPoints+1) = selectedPointEl + maxElevation;
  
if plotProjections
    figure;
    scatter(xTransMercator, yTransMercator);
    hold on,
    scatter(xTransMercator(selectedPointIdx), ...
            yTransMercator(selectedPointIdx), ...
           'r', 'filled', 'LineWidth', 4 );      
    title('Tranverse Mercator projection');
  end    
  

%% Find the points along the elevation

  meridianPath    = [];
  pointNumber     = 0;

  remainingPoints = horzcat(xTransMercator, yTransMercator);
  currentPointIdx = selectedPointIdx;
  cDist = 0;
  
  while validDistance(cDist)
  
    currentPoint    = remainingPoints(currentPointIdx, :);
    remainingPoints(currentPointIdx, :) = inf;
    
    remainingIdx    = remainingPoints(:, 2) > currentPoint(2);
    remainingPoints(~remainingIdx, :) = inf;
    
    pointNumber   = pointNumber+1;
    [cDist, cIdx] = pdist2(remainingPoints, currentPoint, 'euclidean', 'Smallest', 1);
    
    meridianPath(pointNumber).idxPointA = currentPointIdx;
    meridianPath(pointNumber).idxPointB = cIdx;
    meridianPath(pointNumber).distAB    = cDist;
    
    % plot the walk
    if debugMode
      hold on,
      scatter(xTransMercator(cIdx), yTransMercator(cIdx), ...
              'g', 'filled', 'LineWidth', 6 );
      pause();
    end
    
    if beVerbose
      fprintf('Point %d->%d dist %3.2f...\n', currentPointIdx, cIdx, cDist);
    end
    
    % update
    currentPointIdx = cIdx;      
    
  end
  
  % Wipe out duplicated origin
  idxMeridian = [meridianPath.idxPointA];
  meridianPath(idxMeridian == numPoints+1) = [];

  
  % sort the points to close the loop
  idxMeridian       = [meridianPath.idxPointA];
  meridianAzimuth   = currentAzimuth(idxMeridian)+ pi;
  meridianElevation = currentElevation(idxMeridian);
  
  % first half
  idxWithinAz  = find(meridianAzimuth < pi);
  [~, idxA]    = sort(meridianElevation(idxWithinAz), 'ascend');
  idxFirstHalf = idxWithinAz(idxA);

  % second half
  idxWithinAz   = find(meridianAzimuth >= pi);
  [~, idxB]     = sort(meridianElevation(idxWithinAz), 'descend');  
  idxSecondHalf = idxWithinAz(idxB);
 
  idxSortedMeridian  = vertcat(idxFirstHalf, idxSecondHalf);  
  meridianPathSorted = meridianPath;
  
  for meridianIdx = 1:numel(idxSortedMeridian)
    meridianPathSorted(meridianIdx) = meridianPath(idxSortedMeridian(meridianIdx));
  end

  meridianPath = meridianPathSorted;
  
  
  %% Plot the paths if needed
if plotPaths
  
    idxParallel = [parallelPath.idxPointA];
    idxMeridian = [meridianPath.idxPointA];

    figureLeft,
    scatter(currentAzimuth, currentElevation)
    hold on,
    scatter(currentAzimuth(idxMeridian), currentElevation(idxMeridian), 'g', 'filled', 'LineWidth', 4 );
    hold on,
    scatter(currentAzimuth(idxParallel), currentElevation(idxParallel), 'c', 'filled', 'LineWidth', 4 );
    hold on,
    scatter(currentAzimuth(selectedPointIdx), currentElevation(selectedPointIdx), 'r', 'filled', 'LineWidth', 4 );
    hold on,
    plot(currentAzimuth(idxMeridian), currentElevation(idxMeridian), 'LineWidth', 0.5 );
    hold on,
    plot(currentAzimuth(idxParallel), currentElevation(idxParallel), 'LineWidth', 0.5 );
    title('Cyan azimuth');
    
end