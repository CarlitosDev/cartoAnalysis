function  [parallelPath, meridianPath] = getMeridianAndParallelV2(pointsPosition, selectedPointIdx)

% GETMERIDIANANDPARALLELV2 find paths along the azimuth and elevation for a
% given point.
% 
%

%% set up
  % Center points
  doPointCentering = false;
  beVerbose        = true;
  plotPaths        = false;%true;%
  
  validDistance    = @(d) ~(isempty(d) || isnan(d) || isinf(d));  
  
  azimuthRange     = 10*pi/180;
  elevationRange   = 20*pi/180;

%% Center the data and move to spherical coordinates

    if doPointCentering 
        % selected point as origin >> will result in 
        % selectedPointAz =3.1416% pi
        % selectedPointEl =1.5708% pi/2        
        pointsCentroid = pointsCentered(selectedPointIdx, :);
    else
      %  Centroid centering 
      pointsCentroid = mean(pointsPosition, 1);       
    end

    pointsCentered = ...
        bsxfun(@minus, pointsPosition, pointsCentroid);
          
    % Move along the azimuth. Go to spherical coordinates to sort the points by
    % azimuth, so we can move the camera along the horizontal axis.
    [currentAzimuth, currentElevation, currentRadius] = ...
      cart2sph(pointsCentered(:,1), ...
      pointsCentered(:,2), ...
      pointsCentered(:,3));
                                                   
    currentAzimuth   = currentAzimuth   + pi;
    currentElevation = currentElevation + 0.5*pi;

    
%% get the p.o.i and shift the azimuth so the p.o.i gets 0 rads

    selectedPointAz  = currentAzimuth(selectedPointIdx);
    selectedPointEl  = currentElevation(selectedPointIdx);    

    shiftedAz = currentAzimuth-selectedPointAz  ;
    idxShift  = currentAzimuth < selectedPointAz;
    shiftedAz(idxShift) = shiftedAz(idxShift)+(2*pi);

    currentAzimuth   = shiftedAz;

    maxAzimuth     = selectedPointAz + azimuthRange;
    minAzimuth     = selectedPointAz - azimuthRange;

    maxElevation   = selectedPointEl + elevationRange;
    minElevation   = selectedPointEl - elevationRange;
      
    validElevation = currentElevation < maxElevation & ...
                     currentElevation > minElevation;
    

%% azimuth

  parallelPath = [];
  pointNumber  = 0;

  currentPointIdx = selectedPointIdx;
  cDist = 0;
  
  while validDistance(cDist)
  
    currentPoint   = pointsCentered(currentPointIdx, :);
    currentPointAz = currentAzimuth(currentPointIdx);   
    
    validAzimuth   = currentAzimuth > currentPointAz;
    validAzEl      = validElevation & validAzimuth;
    remainingIdx   = find(validAzEl);
    
    pointNumber    = pointNumber+1;
    [cDist, cIdx]  = pdist2(pointsCentered(remainingIdx, :), currentPoint, 'euclidean', 'Smallest', 1);
    
    parallelPath(pointNumber).idxPointA = currentPointIdx;
    parallelPath(pointNumber).idxPointB = remainingIdx(cIdx);
    parallelPath(pointNumber).distAB    = cDist;
    
    if beVerbose
      fprintf('Point %d->%d azimuth %3.2f...\n', currentPointIdx, remainingIdx(cIdx), currentPointAz);
    end
    
    % update
    currentPointIdx = remainingIdx(cIdx);
    
  end
  
%%

meridianPath = [];
return;
    
%% elevation

  meridianPath    = [];
  pointNumber     = 0;

  remainingPoints = horzcat(plateCarreeX, plateCarreeY);
  currentPointIdx = selectedPointIdx;
  cDist = 0;
  
  % half of the orange (STILL MISSING THE LOWER QUADRANT FROM 0, PI/2)
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
    
    if beVerbose
      fprintf('Point %d->%d dist %3.2f...\n', currentPointIdx, cIdx, cDist);
    end
    
    % update
    currentPointIdx = cIdx;
    
  end
  
  % the other half
  % let's load the points again, remove the last pair and shift the last
  % point to the other side by adding pi
    remainingPoints = horzcat(plateCarreeX, plateCarreeY);
    currentPointIdx = meridianPath(pointNumber).idxPointA;
    remainingPoints(currentPointIdx, 1) = remainingPoints(currentPointIdx, 1) + pi;
    pointNumber     = pointNumber - 1; 
    cDist = 0;
  
  while validDistance(cDist)
  
    currentPoint    = remainingPoints(currentPointIdx, :);
    remainingPoints(currentPointIdx, :) = inf;
    
    % move downwards
    remainingIdx    = remainingPoints(:, 2) < currentPoint(2);
    remainingPoints(~remainingIdx, :) = inf;
    
    pointNumber   = pointNumber+1;
    [cDist, cIdx] = pdist2(remainingPoints, currentPoint, 'euclidean', 'Smallest', 1);
    
    meridianPath(pointNumber).idxPointA = currentPointIdx;
    meridianPath(pointNumber).idxPointB = cIdx;
    meridianPath(pointNumber).distAB    = cDist;
    
    if beVerbose
      fprintf('Point %d->%d dist %3.2f...\n', currentPointIdx, cIdx, cDist);
    end
    
    % update
    currentPointIdx = cIdx;
    
  end    
    
  
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