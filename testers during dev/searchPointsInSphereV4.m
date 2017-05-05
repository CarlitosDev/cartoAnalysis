% Generate a unit sphere
    numPoints = 90; % that's per dimension

    maxNumPoints = 2000;
    % let's say we pick a point 
    selectedPointIdx = 83; %randi(numPoints, 1);

    theta=linspace(0,2*pi,numPoints);
    phi=linspace(0,pi,numPoints);
    rho=1;

    [meshTheta,meshPhi] = meshgrid(theta,phi);

    x=rho*sin(meshPhi).*cos(meshTheta);
    y=rho*sin(meshPhi).*sin(meshTheta);
    z=rho*cos(meshPhi);

    % get a unique list of points

    allPoints = unique([x(:), y(:), z(:)], 'rows');
    [numRows, ~] = size(allPoints);

    % randomise the points
    idxPoints = randi(numRows, [maxNumPoints, 1]);

    allPoints = allPoints(idxPoints, :);

        
%% center the data to the selected point

%pointsCentroid = allPoints(selectedPointIdx, :);
pointsCentroid = mean(allPoints, 1);   
pointsCentered = ...
    bsxfun(@minus, allPoints, pointsCentroid);
  
[numPoints, ~] = size(pointsCentered);


%%

figureLeft,
plot3(pointsCentered(:, 1),pointsCentered(:, 2), pointsCentered(:, 3), ...
'o', 'MarkerSize' , 8, ...
'MarkerFaceColor',[.49 1 .63]);

hold on,

plot3(pointsCentered(selectedPointIdx, 1), ...
      pointsCentered(selectedPointIdx, 2), ...
      pointsCentered(selectedPointIdx, 3), ...
      'o', 'MarkerSize' , 20, ...
      'MarkerFaceColor',[1 1 0]);
  






%% go to spherical

[currentAzimuth, currentElevation, currentRadius] = ...
  cart2sph(pointsCentered(:,1), ...
  pointsCentered(:,2), ...
  pointsCentered(:,3));


  %%
  
 [parallelPath, meridianPath] = getMeridianAndParallelMercator(currentAzimuth, currentElevation, selectedPointIdx)
  


%% new devel

idxParallel = [parallelPath.idxPointA];

plot3(pointsCentered(idxParallel, 1), ...
      pointsCentered(idxParallel, 2), ...
      pointsCentered(idxParallel, 3), ...
      'o', 'MarkerSize' , 10, ...
      'MarkerFaceColor',[0 1 1]);

plot3(pointsCentered(idxParallel, 1), ...
      pointsCentered(idxParallel, 2), ...
      pointsCentered(idxParallel, 3), ...
      'MarkerFaceColor',[255,127,80]/255, ...
      'LineWidth', 1.5 );
    
idxMeridian = [meridianPath.idxPointA];

plot3(pointsCentered(idxMeridian, 1), ...
      pointsCentered(idxMeridian, 2), ...
      pointsCentered(idxMeridian, 3), ...
      'o', 'MarkerSize' , 10, ...
      'MarkerFaceColor',[255 140 0]/255);    

plot3(pointsCentered(idxMeridian, 1), ...
      pointsCentered(idxMeridian, 2), ...
      pointsCentered(idxMeridian, 3), ...
      'MarkerFaceColor',[1 0 0], ...
      'LineWidth', 1.5 );
    
%%
% 
% 
% paralX = currentAzimuth;
% paralY = currentElevation;
% 
% figureLeft,
% scatter(paralX, paralY)
% hold on,
% scatter(paralX(idxMeridian), paralY(idxMeridian), 'g', 'filled', 'LineWidth', 4 )
% hold on,
% scatter(paralX(idxParallel), paralY(idxParallel), 'c', 'filled', 'LineWidth', 4 )
% hold on,
% scatter(paralX(selectedPointIdx), paralY(selectedPointIdx), 'r', 'filled', 'LineWidth', 4 )
% hold on,
% plot(paralX(idxMeridian), paralY(idxMeridian), 'LineWidth', 0.5 )
% hold on,
% plot(paralX(idxParallel), paralY(idxParallel), 'LineWidth', 0.5 )