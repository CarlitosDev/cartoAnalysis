function plotSignalsIn3D(hAxTime, timeLength, numSamples, ...
  currentData, pathType, dataFromPath, pointsList, doScaling)

%% inliners
scaleFrom0To1      = @(x) (x-min(x))/(max(x)-min(x));
vectorToStr        = @(x) cellstr(num2str(x));

%%

timeAxes      = linspace(0, timeLength, numSamples);
cumDistPath   = [0, cumsum([dataFromPath.distAB])];
xAxisTickName = vectorToStr(pointsList);

numPathPoints = numel(dataFromPath);

% build a x-axis to place the signal
idx = 1;
currentXAxis = zeros([1, numSamples]);
currentSignal= currentData(idx, :);
if doScaling
  currentTimeValues = scaleFrom0To1(currentSignal);
else
  currentTimeValues = currentSignal;
end

plot3(hAxTime, currentXAxis, ...
  timeAxes, currentTimeValues, 'r', ...
  'LineWidth', 1.2);
hold(hAxTime, 'on');

% plot the remaining signals
for idx = 1:numPathPoints-1
  
  % shift the axis by the distance between points
  currentXAxis = currentXAxis + dataFromPath(idx).distAB;
  
  currentSignal = currentData(idx+1, :);
  if doScaling
    currentTimeValues = scaleFrom0To1(currentSignal);
  else
    currentTimeValues = currentSignal;
  end
  
  plot3(hAxTime , currentXAxis      , ...
    timeAxes, currentTimeValues, 'b', ...
    'LineWidth', 1.0);
  
  hold(hAxTime, 'on');
  
end

hAxTime.View = [-65, 80];
axis(hAxTime, 'tight');
axis(hAxTime, 'ij');

% set point-numbers
hAxTime.XTick      = cumDistPath;
hAxTime.XTickLabel = xAxisTickName;
hAxTime.FontAngle  = 'italic';
hAxTime.FontSize   = 7;
grid(hAxTime, 'on');
parallelTitle = sprintf('%d Points %s', numel(dataFromPath), pathType);
title(hAxTime , parallelTitle)
ylabel(hAxTime, 't(seconds)');