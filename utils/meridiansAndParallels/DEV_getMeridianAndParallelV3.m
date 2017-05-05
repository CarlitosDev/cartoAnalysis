function  [parallelPath, meridianPath] = DEV_getMeridianAndParallelV3(pointsCentered, selectedPointIdx)
% GETMERIDIANANDPARALLELV3 find paths along the azimuth and elevation
% of a (kind of spherical) dataset for a given point.
% 

% The current solution, draws N-points parallel and  meridian starting from
% the selected point and calculates the shortest euclidean distances for
% every point within the lines to the dataset.
%
%
% Carlos Aguilar - 6th November 2k16

%% set up
  
  % Center points
  beVerbose        = true;
  plotPaths        = false;%true;%
  
  validDistance    = @(d) ~(isempty(d) || isnan(d) || isinf(d));  
  
  azimuthRange     = 20*pi/180;
  elevationRange   = 20*pi/180;
  
  % Number of points along the path 
  numPathPoints = 100;    
  
  pointsAlongTheParallel = zeros(1, numPathPoints);
  pointsAlongTheMeridian = zeros(1, numPathPoints);
   
   
    %%
    % Move along the azimuth. Go to spherical coordinates to sort the points by
    % azimuth, so we can move the camera along the horizontal axis.
    [currentAzimuth, currentElevation, currentRadius] = ...
      cart2sph(pointsCentered(:,1), ...
      pointsCentered(:,2), ...
      pointsCentered(:,3));

    selectedPointAz  = currentAzimuth(selectedPointIdx);
    selectedPointEl  = currentElevation(selectedPointIdx);  
   
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
    radiusWithinElevation = currentRadius(validElevation);
    maxRadiusElevation    = 1.2*max(radiusWithinElevation);

    radiusWithinAzimuth = currentRadius(validAzimuth);
    maxRadiusAzimuth    = 1.2*max(radiusWithinAzimuth);

    %%
    % How to build the azimuth?
    fullAzimuth = linspace(0, 2*pi, numPathPoints);

    % Get the parallel in cartesian coords
    [xParallel, yParallel, zParallel] = ...
      sph2cart(fullAzimuth, selectedPointEl, maxRadiusElevation);
    % why z is always one point...?
    zParallel = zParallel*ones(size(xParallel));

figure,
    plot3(pointsCentered(:,1), ...
      pointsCentered(:,2), ...
      pointsCentered(:,3), ...
'o', 'MarkerSize' , 12, ...
'MarkerFaceColor',[.49 1 .63]);
hold on
plot3(xParallel, yParallel, zParallel, 'o')


for idxParallel=1:numPathPoints
  
    currentPoint = [xParallel(idxParallel), ...
        yParallel(idxParallel), ...
        zParallel(idxParallel)];
    
    % TO-DO: Just pass the ones within the elevation range
    [cDist, cIdx]  = pdist2(pointsCentered, currentPoint, 'euclidean', 'Smallest', 1);
    fprintf('current distance...%3.2f\n', cDist);
       
    pointsAlongTheParallel(idxParallel) = cIdx;
    hold on,
    plot3(pointsCentered(cIdx,1), ...
      pointsCentered(cIdx,2), ...
      pointsCentered(cIdx,3), ...
          'o', 'MarkerSize' , 12, ...
        'MarkerFaceColor',[0.1 0.2 0.9]);
    pause(0.1);
  
end

    
%% meridian
numPointsPerSide = ceil(numPathPoints/2);
fullElevation = linspace(-pi/2, pi/2, numPointsPerSide);

%% Part A
[xMeridian, yMeridian, zMeridian] = ...
  sph2cart(selectedPointAz-pi, fullElevation, maxRadiusAzimuth);

for idxMeridian=1:numPointsPerSide
  
    currentPoint = [xMeridian(idxMeridian), ...
        yMeridian(idxMeridian), ...
        zMeridian(idxMeridian)];
    
    % TO-DO: Just pass the ones within the elevation range
    [cDist, cIdx]  = pdist2(pointsCentered, currentPoint, 'euclidean', 'Smallest', 1);
    fprintf('current distance...%3.2f\n', cDist);
       
    pointsAlongTheParallel(idxParallel) = cIdx;
    hold on,
    plot3(pointsCentered(cIdx,1), ...
      pointsCentered(cIdx,2), ...
      pointsCentered(cIdx,3), ...
          'o', 'MarkerSize' , 12, ...
        'MarkerFaceColor',[0.9 0.2 0.2]);
    pause(0.1);
  
end


%% Part B
[xMeridian, yMeridian, zMeridian] = ...
  sph2cart(selectedPointAz, fullElevation, maxRadiusAzimuth);

for idxMeridian=1:numPointsPerSide
  
    currentPoint = [xMeridian(idxMeridian), ...
        yMeridian(idxMeridian), ...
        zMeridian(idxMeridian)];
    
    % TO-DO: Just pass the ones within the elevation range
    [cDist, cIdx]  = pdist2(pointsCentered, currentPoint, 'euclidean', 'Smallest', 1);
    fprintf('current distance...%3.2f\n', cDist);
       
    pointsAlongTheParallel(numPathPoints-idxParallel+1) = cIdx;
    hold on,
    plot3(pointsCentered(cIdx,1), ...
      pointsCentered(cIdx,2), ...
      pointsCentered(cIdx,3), ...
          'o', 'MarkerSize' , 12, ...
        'MarkerFaceColor',[0.9 0.2 0.2]);
    pause(0.1);
  
end


  
%% rearrange the parallel
% drop dups and zeroes

pointsAlongTheParallel = unique(pointsAlongTheParallel, 'stable');
idxToRemove  = pointsAlongTheParallel==0;
pointsAlongTheParallel(idxToRemove) = [];
numParallelPoints = numel(pointsAlongTheParallel);

  parallelPath = [];
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

    parallelPath(pointNumber).idxPointA = idxA;
    parallelPath(pointNumber).idxPointB = idxB;
    parallelPath(pointNumber).distAB    = cDist;
  end
  
  
%% rearrange the meridian
% drop dups and zeroes

pointsAlongTheMeridian = unique(pointsAlongTheMeridian, 'stable');
idxToRemove  = pointsAlongTheMeridian==0;
pointsAlongTheMeridian(idxToRemove) = [];
numMeridianPoints = numel(pointsAlongTheMeridian);

  meridianPath = [];
  for pointNumber=1:numMeridianPoints -1

    idxA = pointsAlongTheParallel(pointNumber);
    idxB = pointsAlongTheParallel(1+pointNumber);

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


  %% Plot the paths if needed
  
if plotPaths

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