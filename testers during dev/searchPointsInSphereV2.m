% Generate a unit sphere
numPoints = 100; % that's per dimension

maxNumPoints = 2000;
% let's say we pick a point 
selectedPointIdx = 13; %randi(numPoints, 1);

theta=linspace(0,2*pi,numPoints);
phi=linspace(0,pi,numPoints);
rho=1;

[meshTheta,meshPhi]=meshgrid(theta,phi);

x=rho*sin(meshPhi).*cos(meshTheta);
y=rho*sin(meshPhi).*sin(meshTheta);
z=rho*cos(meshPhi);

figure
% mesh(x,y,z)
% hold on

%% get a unique list of points

allPoints = unique([x(:), y(:), z(:)], 'rows');
[numRows, ~] = size(allPoints);

% randomise the points
idxPoints = randi(numRows, [maxNumPoints, 1]);

% xAll = allPoints(idxPoints, 1);
% yAll = allPoints(idxPoints, 2);
% zAll = allPoints(idxPoints, 3);

allPoints = allPoints(idxPoints, :);


%% center the data to the selected point

pointsCentroid = mean(allPoints, 1);
pointsCentered = ...
    bsxfun(@minus, allPoints, pointsCentroid);
  
[numPoints, ~] = size(pointsCentered);

currentPointXYZ  = pointsCentered(selectedPointIdx, :); 

%%


%% let's just do the translation
% use bsxfun(@minus, pointsPosition, pointsCentroid);
for idx=1:3
  pointsCentered(:, idx) = pointsCentered(:, idx)  - currentPointXYZ(idx);
end

% by doing the translation, when moving to sphericals:
% selectedPointAz =3.1416% pi
% selectedPointEl =1.5708% pi/2

% by rotating...
% selectedPointAz =5.5885 ??? % pi >> 1.3656
% selectedPointEl =1.5708% pi/2

% phi -> x-axis
% psi -> y-axis

% phi = -pi;
% psi = -pi/2;
% chi = -pi;

% a = [cos(psi)*cos(chi), -cos(psi)*sin(chi), sin(psi); ...
%   
% cos(phi)*sin(chi)+sin(phi)*sin(psi)*cos(chi), ...
% cos(phi)*cos(chi)-sin(phi)*sin(psi)*sin(chi), ...
% -sin(phi)*cos(psi); ...
% 
% sin(phi)*sin(chi)-cos(phi)*sin(psi)*cos(chi), ...
% sin(phi)*cos(chi)+cos(phi)*sin(psi)*sin(chi), ...
% cos(phi)*cos(psi)];
% 
% % Xrt = a(1,1)*x+a(1,2)*y+a(1,3)*z+Lx;
% idxTransform = 1;
% pointsCentered(:, idxTransform) = ...
%   a(idxTransform, 1).*pointsCentered(:, 1) + ...
%   a(idxTransform, 2).*pointsCentered(:, 2) + ...
%   a(idxTransform, 3).*pointsCentered(:, 3);
%   
% % Yrt = a(2,1)*x+a(2,2)*y+a(2,3)*z+Ly;
% idxTransform = 2;
% pointsCentered(:, idxTransform) = ...
%   a(idxTransform, 1).*pointsCentered(:, 1) + ...
%   a(idxTransform, 2).*pointsCentered(:, 2) + ...
%   a(idxTransform, 3).*pointsCentered(:, 3);
% 
% % Zrt = a(3,1)*x+a(3,2)*y+a(3,3)*z+Lz;
% idxTransform = 3;
% pointsCentered(:, idxTransform) = ...
%   a(idxTransform, 1).*pointsCentered(:, 1) + ...
%   a(idxTransform, 2).*pointsCentered(:, 2) + ...
%   a(idxTransform, 3).*pointsCentered(:, 3);



%a = ones(numPoints, 3)*([1,2,3]'*ones(1,3))'/3;

% translationMatrix = -1*ones(numPoints, 3);
% 
% idxTranslation = 1;
% translationMatrix(:, idxTranslation) = ...
%   translationMatrix(:, idxTranslation)*currentPointXYZ(idxTranslation);
% 
% idxTranslation = 1 + idxTranslation;
% translationMatrix(:, idxTranslation) = ...
%   translationMatrix(:, idxTranslation)*currentPointXYZ(idxTranslation);
% 
% idxTranslation = 1 + idxTranslation;
% translationMatrix(:, idxTranslation) = ...
%   translationMatrix(:, idxTranslation)*currentPointXYZ(idxTranslation);
% 
%  
% pointsCentered = pointsCentered + translationMatrix;

%% go to spherical

[currentAzimuth, currentElevation, currentRadius] = ...
  cart2sph(pointsCentered(:,1), ...
  pointsCentered(:,2), ...
  pointsCentered(:,3));

currentAzimuth   = currentAzimuth   + pi;
currentElevation = currentElevation + 0.5*pi;

%% rotate and translate the matrix to the current point (set as origin)
  



% selectedPointAz  = currentAzimuth(selectedPointIdx);
% selectedPointEl  = currentElevation(selectedPointIdx);
% 
% currentPointXYZ  = pointsCentered(selectedPointIdx, :);

% M = makehgtform('translate', -1*currentPointXYZ);

% psi = 0;
% chi = 0;
% phi = 0;
% x   = currentPointXYZ(1);
% y   = currentPointXYZ(2);
% z   = currentPointXYZ(3);
% Lx  = -1*currentPointXYZ(1);
% Ly  = -1*currentPointXYZ(2);
% Lz  = -1*currentPointXYZ(3);
% 
% [Xrt, Yrt, Zrt] = transformEulerAngles(psi, chi, phi, Lx, Ly, Lz, x, y, z)

% that's the function transformEulerAngles
%
% a = [cos(psi)*cos(chi), -cos(psi)*sin(chi), sin(psi); ...
%   
% cos(phi)*sin(chi)+sin(phi)*sin(psi)*cos(chi), ...
% cos(phi)*cos(chi)-sin(phi)*sin(psi)*sin(chi), ...
% -sin(phi)*cos(psi); ...
% 
% sin(phi)*sin(chi)-cos(phi)*sin(psi)*cos(chi), ...
% sin(phi)*cos(chi)+cos(phi)*sin(psi)*sin(chi), ...
% cos(phi)*cos(psi)];
% 
% 
% Xrt = a(1,1)*x+a(1,2)*y+a(1,3)*z+Lx;
% Yrt = a(2,1)*x+a(2,2)*y+a(2,3)*z+Ly;
% Zrt = a(3,1)*x+a(3,2)*y+a(3,3)*z+Lz;

%%




%% 

plot3(pointsCentered(:, 1),pointsCentered(:, 2), pointsCentered(:, 3), ...
'o', 'MarkerSize' , 12, ...
'MarkerFaceColor',[.49 1 .63]);

hold on,

plot3(pointsCentered(selectedPointIdx, 1), ...
      pointsCentered(selectedPointIdx, 2), ...
      pointsCentered(selectedPointIdx, 3), ...
      'o', 'MarkerSize' , 20, ...
      'MarkerFaceColor',[1 1 0]);


%% new devel

% equirectangular projection for the parallels
[parallelPath, meridianPath] = ...
  getMeridianAndParallel(currentAzimuth, currentElevation, selectedPointIdx);

idxParallel = [parallelPath.idxPointA];

plot3(pointsCentered(idxParallel, 1), ...
      pointsCentered(idxParallel, 2), ...
      pointsCentered(idxParallel, 3), ...
      'o', 'MarkerSize' , 10, ...
      'MarkerFaceColor',[0 1 1]);

plot3(pointsCentered(idxParallel, 1), ...
      pointsCentered(idxParallel, 2), ...
      pointsCentered(idxParallel, 3), ...
      'MarkerFaceColor',[1 0 0], ...
      'LineWidth', 1.5 );
    
idxMeridian = [meridianPath.idxPointA];

plot3(pointsCentered(idxMeridian, 1), ...
      pointsCentered(idxMeridian, 2), ...
      pointsCentered(idxMeridian, 3), ...
      'o', 'MarkerSize' , 10, ...
      'MarkerFaceColor',[1 0 0]);    

plot3(pointsCentered(idxMeridian, 1), ...
      pointsCentered(idxMeridian, 2), ...
      pointsCentered(idxMeridian, 3), ...
      'MarkerFaceColor',[1 0 0], ...
      'LineWidth', 1.5 );
    
%%
% 
% 
paralX = currentAzimuth;
paralY = currentElevation;

figureLeft,
scatter(paralX, paralY)
hold on,
scatter(paralX(idxMeridian), paralY(idxMeridian), 'g', 'filled', 'LineWidth', 4 )
hold on,
scatter(paralX(idxParallel), paralY(idxParallel), 'c', 'filled', 'LineWidth', 4 )
hold on,
scatter(paralX(selectedPointIdx), paralY(selectedPointIdx), 'r', 'filled', 'LineWidth', 4 )
hold on,
plot(paralX(idxMeridian), paralY(idxMeridian), 'LineWidth', 0.5 )
hold on,
plot(paralX(idxParallel), paralY(idxParallel), 'LineWidth', 0.5 )