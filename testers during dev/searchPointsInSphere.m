% Generate a unit sphere
numPoints = 60; % that's per dimension

theta=linspace(0,2*pi,numPoints);
phi=linspace(0,pi,numPoints);
rho=1;

[meshTheta,meshPhi]=meshgrid(theta,phi);


x=rho*sin(meshPhi).*cos(meshTheta);
y=rho*sin(meshPhi).*sin(meshTheta);
z=rho*cos(meshPhi);

figure
mesh(x,y,z)
hold on

%% get a unique list of points

xAll = x(:);
yAll = y(:);
zAll = z(:);

allPoints = unique([x(:), y(:), z(:)], 'rows');
[numRows, ~] = size(allPoints);

% randomise the points
maxNumPoints = 300;
idxPoints = randi(numRows, [maxNumPoints, 1]);


xAll = allPoints(idxPoints, 1);
yAll = allPoints(idxPoints, 2);
zAll = allPoints(idxPoints, 3);

plot3(xAll,yAll,zAll, ...
'o', 'MarkerSize' , 12, ...
'MarkerFaceColor',[.49 1 .63]);


%%

pointsCentered = allPoints(idxPoints, :);

[currentAzimuth, currentElevation, currentRadius] = ...
  cart2sph(pointsCentered(:,1), ...
  pointsCentered(:,2), ...
  pointsCentered(:,3));

currentAzimuth   = currentAzimuth   + pi;
currentElevation = currentElevation + 0.5*pi;

% let's say we pick a point 
selectedPointIdx = 153; %randi(numPoints, 1);
selectedPointAz  = currentAzimuth(selectedPointIdx);
selectedPointEl  = currentElevation(selectedPointIdx);


%%

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


%% create meridian rings

[xTemp, yTemp, zTemp] = sph2cart(minAzimuth-pi, fullElevation, maxRadius);
meridianLeftA  = [xTemp; yTemp; zTemp];
[xTemp, yTemp, zTemp] = sph2cart(minAzimuth   , fullElevation, maxRadius);
meridianLeftB  = [xTemp; yTemp; zTemp];
meridianLeft   = [meridianLeftA, meridianLeftB];

% dark side of the moon...
% if maxAzimuth >= pi
% end
    
[xTemp, yTemp, zTemp] = sph2cart(maxAzimuth-pi, fullElevation, maxRadius);
meridianRightA = [xTemp; yTemp; zTemp];
[xTemp, yTemp, zTemp] = sph2cart(maxAzimuth   , fullElevation, maxRadius);
meridianRightB = [xTemp; yTemp; zTemp];
meridianRight  = [meridianRightA, meridianRightB];

hold on,
plot3(meridianLeft(1,:), meridianLeft(2,:), meridianLeft(3,:));

hold on,
plot3(meridianRight(1,:), meridianRight(2,:), meridianRight(3,:));


%% create parallel rings

% x = r .* cos(elevation) .* cos(azimuth)
% y = r .* cos(elevation) .* sin(azimuth)
% z = r .* sin(elevation)

% that will make 'z' as a single number...

[xTemp, yTemp, zTemp] = sph2cart(fullAzimuth, minElevation-0.5*pi, maxRadius);
zTemp = zTemp.*ones(size(yTemp));
parallelDown = [xTemp; yTemp; zTemp];



[xTemp, yTemp, zTemp] = sph2cart(fullAzimuth, maxElevation-0.5*pi, maxRadius);
zTemp = zTemp.*ones(size(yTemp));
parallelUp = [xTemp; yTemp; zTemp];


hold on,
plot3(parallelDown(1,:), parallelDown(2,:), parallelDown(3,:));

hold on,
plot3(parallelUp(1,:), parallelUp(2,:), parallelUp(3,:));

%% 
hold on,

plot3(allPoints(idxSorted, 1), ...
      allPoints(idxSorted, 2), ...
      allPoints(idxSorted, 3), ...
      'o', 'MarkerSize' , 12, ...
      'MarkerFaceColor',[1 0 0]);

hold on,

plot3(pointsCentered(selectedPointIdx, 1), ...
      pointsCentered(selectedPointIdx, 2), ...
      pointsCentered(selectedPointIdx, 3), ...
      'o', 'MarkerSize' , 20, ...
      'MarkerFaceColor',[1 1 0]);




%% new devel

% equirectangular projection for the parallels
[parallelPath, meridianPath] = getMeridianAndParallel(currentAzimuth, currentElevation, selectedPointIdx);

idxParallel = [parallelPath.idxPointA];


plot3(pointsCentered(idxParallel, 1), ...
      pointsCentered(idxParallel, 2), ...
      pointsCentered(idxParallel, 3), ...
      'o', 'MarkerSize' , 10, ...
      'MarkerFaceColor',[0 1 1]);
    
    
idxMeridian = [meridianPath.idxPointA];

plot3(pointsCentered(idxMeridian, 1), ...
      pointsCentered(idxMeridian, 2), ...
      pointsCentered(idxMeridian, 3), ...
      'o', 'MarkerSize' , 10, ...
      'MarkerFaceColor',[0 1 1]);    

%%


figureLeft,
scatter(paralX, paralY)
hold on,
scatter(paralX(selectedPointIdx), paralY(selectedPointIdx), 'r', 'filled', 'LineWidth', 4 )



figureRight,
scatter(shiftedAz, shiftedEl)
hold on,
scatter(shiftedAz(selectedPointIdx), shiftedEl(selectedPointIdx), 'r', 'filled', 'LineWidth', 4 )