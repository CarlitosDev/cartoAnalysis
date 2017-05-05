function plotDataAnalysisParallelAndMeridianV4(cathererTable, ...
  doPlotMesh, meshData, ...
  selectedPointIdx, textPosition, parallelPath, meridianPath, pathInfo)

%% PLOTDATAANALYSISPARALLELANDMERIDIANV4
% Similar to PLOTDATAANALYSISPARALLELANDMERIDIANV3 but it does NOT plot the
% mesh over the carto points
%
% Carlos Aguilar  17/04/2017


%% hard-coded parameters
numSamples        = size(cathererTable{1, 'signal'}, 2);
samplingFrequency = 1e3;
timeLength        = numSamples/samplingFrequency;
fCutOff           = 80;

doScaling         = true;

% define colours
meshColour         = [232,255,255]/255;
markersColour      = [1 0 0];
parallelLineColour = [255,127,80] /255;
meridianLineColour = [255 140 0]  /255;
allPointsColour    = [240, 240, 240]/255;


xData = cathererTable{:, 'x'};
yData = cathererTable{:, 'y'};
zData = cathererTable{:, 'z'};

%% inliners

scaleFrom0To1      = @(x) (x-min(x))/(max(x)-min(x));
vectorToStr        = @(x) cellstr(num2str(x));
%% Input data

% bundle all the points to plot them as small spheres and add the point
% number
idxAllPoints      = unique([pathInfo.allParallelPointsIdx , ...
  pathInfo.allMeridianPointsIdx]);
numSelectedPoints = numel(idxAllPoints);


%% Plot mesh and selected points

h = figure;
set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);

subplot(2,3,[4,1]);

if doPlotMesh
  %     the new mesh is so dense that I can't really see the points. Switch it
  %     off for the time being.
  trimesh(meshData.faces, ...
    meshData.X, meshData.Y, meshData.Z, ...
    'FaceColor', 'none', ...
    'EdgeColor', meshColour); %#ok<UNRCH>
else
  fprintf('\n...Switching off mesh viz\n');
end

% Add points Id
for i=1:numSelectedPoints
  hold on;
  text(textPosition(idxAllPoints(i), 1), ...
    textPosition(idxAllPoints(i), 2), ...
    textPosition(idxAllPoints(i), 3), ...
    sprintf('%d', cathererTable.pointsId(idxAllPoints(i))));
end


set(0, 'DefaultTextInterpreter', 'none');

xLabelString = ...
  sprintf('Point selected %d', ...
  cathererTable.pointsId(selectedPointIdx));
xlabel(xLabelString);

%% plot the 3D mesh, the parallels and meridians and the selected points

% plot all the points belonging to the parallels and meridians
plot3(xData, yData, zData, ...
  'o', 'MarkerSize' , 3, ...
  'MarkerFaceColor' , allPointsColour, ...
  'MarkerEdgeColor' , 'none');

% plot all the points belonging to the parallels and meridians
plot3(xData(idxAllPoints), ...
  yData(idxAllPoints), ...
  zData(idxAllPoints), ...
  'o', 'MarkerSize'  , 10, ...
  'MarkerFaceColor'  , markersColour);

% plot the selected point in yellow colour
plot3(xData(selectedPointIdx), ...
  yData(selectedPointIdx), ...
  zData(selectedPointIdx), ...
  'o', 'MarkerSize'      , 20, ...
  'MarkerFaceColor'      ,[1 1 0]);

% plot meridian points
plot3(xData(pathInfo.allMeridianPointsIdx) , ...
  yData(pathInfo.allMeridianPointsIdx) , ...
  zData(pathInfo.allMeridianPointsIdx) , ...
  'MarkerFaceColor', meridianLineColour, ...
  'LineWidth', 1.8 );

% plot parallel points
plot3(xData(pathInfo.allParallelPointsIdx) , ...
  yData(pathInfo.allParallelPointsIdx) , ...
  zData(pathInfo.allParallelPointsIdx) , ...
  'MarkerFaceColor', parallelLineColour, ...
  'LineWidth', 1.8 );

% point the plot to the selected point
cAxes = gca();
axis(cAxes, 'tight')

drawnow;

%% Get the point Id's

% parallel
parallelPointsList    = ...
  cathererTable.pointsId(pathInfo.allParallelPointsIdx);
xAxisParallelTickName = vectorToStr(parallelPointsList);

% meridian
meridianPointsList    = ...
  cathererTable.pointsId(pathInfo.allMeridianPointsIdx);
xAxisMeridianTickName = vectorToStr(meridianPointsList);


%% (2.1) Time-Plot EGM parallel

% get the closed path
idxClosedPath = ...
  [parallelPath(1).idxPointA, parallelPath.idxPointB];

% get the signals bounded to the path
currentData = cathererTable{idxClosedPath, 'signal'};

hAxTime           = subplot(2,3,[2,2]);
numSamples        = size(currentData, 2);
samplingFrequency = 1e3;
timeLength        = numSamples/samplingFrequency;

doScaling  = true;
pathType   = 'parallel';
pointsList = ...
  cathererTable.pointsId(pathInfo.allParallelPointsIdx);

plotSignalsIn3D(hAxTime, timeLength, numSamples, ...
  currentData, pathType, parallelPath, pointsList, doScaling);

%% (2.2)  Plot frequency EGM parallel

hAx = subplot(2,3,[5,5]);
cumDistParallel = [0, cumsum([parallelPath.distAB])];

% first plot the selected point
idx = 1;
idxCurrentParallel = parallelPath(idx).idxPointA;

currentCartoData   = cathererTable{idxCurrentParallel, 'signal'};

[inputPSD, currentFreqAxes] = ...
  getPowerSpectralDensity(currentCartoData, samplingFrequency, ...
  'windowType', 'rectangular', 'windowSize', numSamples, ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', fCutOff);

currentPSDLength   = numel(currentFreqAxes);
currentXAxis       = zeros([1, currentPSDLength]);

if doScaling
  currentFreqValues = scaleFrom0To1(inputPSD);
else
  currentFreqValues = inputPSD; %#ok<UNRCH>
end
plot3(hAx, currentXAxis, ...
  currentFreqAxes, currentFreqValues, ...
  'r', 'LineWidth', 1.0);
hold(hAx, 'on');

for idx = 1:pathInfo.numParallelPoints-1
  
  idxCurrentParallel = parallelPath(idx).idxPointB;
  
  currentCartoData  = cathererTable{idxCurrentParallel, 'signal'};
  
  [inputPSD, currentFreqAxes] = ...
    getPowerSpectralDensity(currentCartoData, samplingFrequency, ...
    'windowType', 'rectangular', 'windowSize', numSamples, ...
    'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', fCutOff);
  
  if doScaling
    currentFreqValues = scaleFrom0To1(inputPSD);
  else
    currentFreqValues = inputPSD; %#ok<UNRCH>
  end
  
  currentPSDLength = numel(currentFreqAxes);
  currentXAxis     = cumDistParallel(idx+1)*ones([1, currentPSDLength]);
  
  plot3(hAx, currentXAxis, ...
    currentFreqAxes, currentFreqValues, ...
    'b', 'LineWidth', 1.0);
  
  hold(hAx, 'on');
  
end

axis(hAx, 'tight')
axis(hAx, 'ij');
hAx.View       = [-65, 80];
hAx.XTick      = cumDistParallel;
hAx.XTickLabel = xAxisParallelTickName;
hAx.FontAngle  = 'italic';
hAx.FontSize   = 8;
grid(hAx, 'on');
ylabel(hAx, 'f(Hertz)');


%% (3.1)  Time-Plot EGM meridian

% get the closed path
idxClosedPath = ...
  [meridianPath(1).idxPointA, meridianPath.idxPointB];

% get the signals bounded to the path
currentData = cathererTable{idxClosedPath, 'signal'};

hAxMrdn           = subplot(2,3,[3,3]);
numSamples        = size(currentData, 2);
samplingFrequency = 1e3;
timeLength        = numSamples/samplingFrequency;

doScaling  = true;
pathType   = 'meridian';
pointsList = ...
  cathererTable.pointsId(pathInfo.allMeridianPointsIdx);

plotSignalsIn3D(hAxMrdn, timeLength, numSamples, ...
  currentData, pathType, meridianPath, pointsList, doScaling);


%% (3.2) Plot frequency EGM parallel
hAxMrdnFreq = subplot(2,3,[6,6]);

cumDistMeridian = [0, cumsum([meridianPath.distAB])];

% first plot the selected point
idx = 1;
idxCurrentMeridian = meridianPath(idx).idxPointA;

currentCartoData   = cathererTable{idxCurrentMeridian, 'signal'};

[inputPSD, currentFreqAxes] = ...
  getPowerSpectralDensity(currentCartoData, samplingFrequency, ...
  'windowType', 'rectangular', 'windowSize', numSamples, ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', fCutOff);

currentPSDLength   = numel(currentFreqAxes);
currentXAxis       = zeros([1, currentPSDLength]);

if doScaling
  currentFreqValues = scaleFrom0To1(inputPSD);
else
  currentFreqValues = inputPSD; %#ok<UNRCH>
end
plot3(hAxMrdnFreq, currentXAxis, ...
  currentFreqAxes, currentFreqValues, ...
  'r', 'LineWidth', 1.0);

hold(hAxMrdnFreq, 'on');

for idx = 1:pathInfo.numMeridianPoints-1
  
  idxCurrentMeridian = meridianPath(idx).idxPointB;
  
  currentCartoData   = cathererTable{idxCurrentMeridian, 'signal'};
  [inputPSD, currentFreqAxes] = ...
    getPowerSpectralDensity(currentCartoData, samplingFrequency, ...
    'windowType', 'rectangular', 'windowSize', numSamples, ...
    'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', fCutOff);
  if doScaling
    currentFreqValues = scaleFrom0To1(inputPSD);
  else
    currentFreqValues = inputPSD; %#ok<UNRCH>
  end
  
  currentPSDLength = numel(currentFreqAxes);
  currentXAxis     = cumDistMeridian(idx+1)*ones([1, currentPSDLength]);
  
  plot3(hAxMrdnFreq, currentXAxis, ...
    currentFreqAxes, currentFreqValues, ...
    'b', 'LineWidth', 1.0);
  
  hold(hAxMrdnFreq, 'on');
  
end

axis(hAxMrdnFreq, 'tight')
axis(hAxMrdnFreq, 'ij');
hAxMrdnFreq.View       = [-65, 80];
hAxMrdnFreq.XTick      = cumDistMeridian;
hAxMrdnFreq.XTickLabel = xAxisMeridianTickName;
hAxMrdnFreq.FontAngle  = 'italic';
hAxMrdnFreq.FontSize   = 8;
grid(hAxMrdnFreq, 'on');
ylabel(hAxMrdnFreq, 'f(Hertz)');

drawnow;