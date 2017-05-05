function  [parallelPath, meridianPath, pathInfo] = getMeridianAndParallelV4(currentAzimuth, currentElevation, currentRadius, selectedPointIdx, pointsCentered)
% GETMERIDIANANDPARALLELV4 find paths along the azimuth and elevation
% of a (kind of spherical) dataset for a given point.
% 

% The current solution, draws N-points parallel and  meridian starting from
% the selected point and calculates the shortest euclidean distances for
% every point within the lines to the dataset.
%
%
% Carlos Aguilar - 8th November 2k16
      pathInfo     = [];
      parallelPath = [];
      meridianPath = [];
      

%% set up
  azimuthRange     = 20*pi/180;
  elevationRange   = 20*pi/180;
  
  % Number of points along the path 
  numPathPoints = 150;
  
  pointsAlongTheParallel = zeros(1, numPathPoints);
  pointsAlongTheMeridian = zeros(1, numPathPoints);
  
  selectedPointAz  = currentAzimuth(selectedPointIdx);
  selectedPointEl  = currentElevation(selectedPointIdx);
  
  % draw the paths along an average radius or use the selected point radius
  % otherwise.
  useAverageRadius = false;
   
%% parallel

% Get the p.o.i and draw a parallel with the same azimuth and elevation
% that goes from 0 to 2*pi.
% Find the closest point in the dataset for each point within the parallel.

    maxAzimuth     = selectedPointAz + azimuthRange;
    minAzimuth     = selectedPointAz - azimuthRange;

    maxElevation   = selectedPointEl + elevationRange;
    minElevation   = selectedPointEl - elevationRange;

    validElevation = currentElevation < maxElevation & ...
      currentElevation > minElevation;

    validAzimuth   = currentAzimuth < maxAzimuth & ...
      currentAzimuth >  minAzimuth;

    % Can we get the radius based on the data?
    if useAverageRadius
      radiusWithinElevation = currentRadius(validElevation); %#ok<UNRCH>
      radiusElevation       = 1.2*max(radiusWithinElevation);
      radiusWithinAzimuth   = currentRadius(validAzimuth);
      radiusAzimuth         = 1.2*max(radiusWithinAzimuth);
    else
      radiusElevation       = currentRadius(selectedPointIdx);
      radiusAzimuth         = currentRadius(selectedPointIdx);      
    end
    
    % generate the parallel
    fullAzimuth = linspace(0, 2*pi, numPathPoints);

    % Get the parallel in cartesian coords
    [xParallel, yParallel, zParallel] = ...
      sph2cart(fullAzimuth, selectedPointEl, radiusElevation);
    % why z is always one point...?
    zParallel = zParallel*ones(size(xParallel));
    
    % get the parallel
    for idxParallel=1:numPathPoints
        
        currentPoint = [xParallel(idxParallel), ...
            yParallel(idxParallel), ...
            zParallel(idxParallel)];
        
        % TO-DO: Just pass the ones within the elevation range
        [~, cIdx]  = pdist2(pointsCentered, currentPoint, ...
                        'euclidean', 'Smallest', 1);        
        
        pointsAlongTheParallel(idxParallel) = cIdx;
        
    end

    
%% meridian

    % Split the calculation into two parts
    numPointsPerSide = ceil(numPathPoints/2);
    fullElevation = linspace(-pi/2, pi/2, numPointsPerSide);

    % Part A
    [xMeridian, yMeridian, zMeridian] = ...
        sph2cart(selectedPointAz-pi, fullElevation, radiusAzimuth);

    for idxMeridian=1:numPointsPerSide

        currentPoint = [xMeridian(idxMeridian), ...
            yMeridian(idxMeridian), ...
            zMeridian(idxMeridian)];

        % TO-DO: Just pass the ones within the elevation range
        [~, cIdx]  = pdist2(pointsCentered, currentPoint, ...
                     'euclidean', 'Smallest', 1);
        pointsAlongTheMeridian(idxMeridian) = cIdx;

    end

    % Part B
    [xMeridian, yMeridian, zMeridian] = ...
      sph2cart(selectedPointAz, fullElevation, radiusAzimuth);

    for idxMeridian=1:numPointsPerSide

        currentPoint = [xMeridian(idxMeridian), ...
            yMeridian(idxMeridian), ...
            zMeridian(idxMeridian)];

        % TO-DO: Just pass the ones within the elevation range
        [~, cIdx]  = pdist2(pointsCentered, currentPoint, ...
                     'euclidean', 'Smallest', 1);

        pointsAlongTheMeridian(numPathPoints-idxMeridian+1) = cIdx;

    end        
    
%% rearrange the parallel    
    % drop dups and zeroes and calculate the distance between pathpoints
    
    pointsAlongTheParallel = unique(pointsAlongTheParallel, 'stable');
    idxToRemove  = pointsAlongTheParallel==0;
    pointsAlongTheParallel(idxToRemove) = [];
        
    % if is empty the algorith is a bit crappy...
    idxSelected = find(pointsAlongTheParallel == selectedPointIdx);        
    
    %rewind the path so the first and last one in the path is the one selected
    pointsAlongTheParallel = [pointsAlongTheParallel(idxSelected:end), ...
                              pointsAlongTheParallel(1:idxSelected-1), ....
                              pointsAlongTheParallel(idxSelected)];    
    
    % save the unique list
    pathInfo.allParallelPointsIdx = pointsAlongTheParallel;
    
    numParallelPoints = numel(pointsAlongTheParallel);
    pathInfo.numParallelPoints = numParallelPoints;

      for pointNumber=1:numParallelPoints-1

        idxA = pointsAlongTheParallel(pointNumber);
        idxB = pointsAlongTheParallel(1+pointNumber);

        pointA = [pointsCentered(idxA, 1), ...
          pointsCentered(idxA, 2), ...
          pointsCentered(idxA, 3)];

        pointB = [pointsCentered(idxB, 1), ...
          pointsCentered(idxB, 2), ...
          pointsCentered(idxB, 3)];

        cDist = pdist2(pointA , pointB, 'euclidean');

        parallelPath(pointNumber).idxPointA = idxA; %#ok<*AGROW>
        parallelPath(pointNumber).idxPointB = idxB;
        parallelPath(pointNumber).distAB    = cDist;
      end      
  
%% rearrange the meridian
    % drop dups and zeroes and calculate the distance between pathpoints

    pointsAlongTheMeridian = unique(pointsAlongTheMeridian, 'stable');
    idxToRemove = pointsAlongTheMeridian==0;
    pointsAlongTheMeridian(idxToRemove) = [];
    
    % if is empty the algorithm is a bit crappy...
    idxSelected = find(pointsAlongTheMeridian == selectedPointIdx);        
    
    %rewind the path so the first and last one in the path is the one selected
    pointsAlongTheMeridian = [pointsAlongTheMeridian(idxSelected:end), ...
                              pointsAlongTheMeridian(1:idxSelected-1), ....
                              pointsAlongTheMeridian(idxSelected)];    
    
    % save the unique list
    pathInfo.allMeridianPointsIdx = pointsAlongTheMeridian;      
    numMeridianPoints = numel(pointsAlongTheMeridian);
    pathInfo.numMeridianPoints = numMeridianPoints;
    
    for pointNumber=1:numMeridianPoints -1
        
        idxA = pointsAlongTheMeridian(pointNumber);
        idxB = pointsAlongTheMeridian(1+pointNumber);
        
        pointA = [pointsCentered(idxA, 1), ...
                  pointsCentered(idxA, 2), ...
                  pointsCentered(idxA, 3)];
        
        pointB = [pointsCentered(idxB, 1), ...
                  pointsCentered(idxB, 2), ...
                  pointsCentered(idxB, 3)];
        
        cDist = pdist2(pointA , pointB, 'euclidean');
        
        meridianPath(pointNumber).idxPointA = idxA;
        meridianPath(pointNumber).idxPointB = idxB;
        meridianPath(pointNumber).distAB    = cDist;
    end  