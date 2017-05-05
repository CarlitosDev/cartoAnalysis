% Generate a unit sphere
numPoints = 80; % that's per dimension

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
selectedPointIdx = 53; %randi(numPoints, 1);

%%

[parallelPath, meridianPath] = getMeridianAndParallelV3(pointsCentered, selectedPointIdx);