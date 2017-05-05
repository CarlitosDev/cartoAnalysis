function plotCartoDataAnalysisSliced(idxSorted, cartoPointAnalysis, meshData, selectedPointsData)


%% hard-coded parameters (review this...)

numSamples        = 2500;
samplingFrequency = 1e3;
timeLength        = numSamples/samplingFrequency;

doScaling         = true;

%% inliners

scaleFrom0To1      = @(x) (x-min(x))/(max(x)-min(x));
distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));

%% Plot mesh and selected points
  selectedPointIdx  = selectedPointsData.selectedIdx;
  numSelectedPoints = numel(selectedPointsData.xData);

  h = figureRight;
  subplot(2,2,[3,1])

  % plot carto points
  plot3(selectedPointsData.xData, ...
        selectedPointsData.yData, ...
        selectedPointsData.zData, ...
        'o', 'MarkerSize' , 12  , ...
        'MarkerFaceColor',[.49 1 .63]);
      
  % selected one in red colour
  hold on;
  plot3(selectedPointsData.xData(selectedPointIdx), ...
        selectedPointsData.yData(selectedPointIdx), ...
        selectedPointsData.zData(selectedPointIdx), ...
        'o', 'MarkerSize' , 20  , ...
        'MarkerFaceColor',[1 0 0]);
  
	% plot mesh
  hold on;
  trimesh(meshData.faces, ...
          meshData.X, meshData.Y, meshData.Z, ...
          'FaceColor', 'none', ...
          'EdgeColor', [0 0 1]);
  
  % Add points Id  
  for i=1:numSelectedPoints
    hold on;
    text(selectedPointsData.textPosition(i, 1), ...
         selectedPointsData.textPosition(i, 2), ...
         selectedPointsData.textPosition(i, 3), ...
         sprintf('%d', selectedPointsData.pointsId(i)));
  end
  
  xlabel(sprintf('%d points selected', numSelectedPoints));
  
  

%% Get DISTANCES between points. Use Euclidean for the time being - use ortoleches later on
  
  currentPoint  = [selectedPointsData.xData(selectedPointIdx), ...
                   selectedPointsData.yData(selectedPointIdx), ...
                   selectedPointsData.zData(selectedPointIdx)];

  allPoints     = [selectedPointsData.xData, ...
                   selectedPointsData.yData, ...
                   selectedPointsData.zData];
                 
  xAxisPosition = distPoint2PointSet(allPoints, currentPoint);
  xAxisTickName = cellstr(num2str(selectedPointsData.pointsId')); %#ok<NASGU>


%% Time-Plot EGM

  hAxTime  = subplot(2,2,[2,2]);
  timeAxes = linspace(0, timeLength, numSamples);

  for idx = 1:numel(idxSorted)

    currentCartoData  = cartoPointAnalysis{idxSorted(idx)};
    if doScaling
        currentTimeValues = scaleFrom0To1(currentCartoData.dataForICA);
    else
        currentTimeValues = currentCartoData.dataForICA; %#ok<UNRCH>
    end

    currentXAxis = repmat(xAxisPosition(idx), [1, numSamples]);

    if idx == selectedPointIdx
    plot3(hAxTime, currentXAxis, ...
          timeAxes, currentTimeValues, 'r', ...
          'LineWidth', 1.2);
    else
    plot3(hAxTime, currentXAxis, ...
          timeAxes, currentTimeValues, 'b', ...
          'LineWidth', 1.0);
    end
    hold(hAxTime, 'on');

  end
  
  hAxTime.View = [-75.26, 62.64];
  
  
  
  % set point-numbers
  [sortedXAxisPosition, idxXAxis] = sort(xAxisPosition);
  hAxTime.XTick      = sortedXAxisPosition;
  hAxTime.XTickLabel = xAxisTickName(idxXAxis);

%% Plot frequency
hAx = subplot(2,2,[4,4]);

for idx = 1:numel(idxSorted)
  
  currentCartoData  = cartoPointAnalysis{idxSorted(idx)};
  currentFreqAxes   = currentCartoData.inputPSDOmega;
  if doScaling
      currentFreqValues = scaleFrom0To1(currentCartoData.inputPSD);
  else
      currentFreqValues = currentCartoData.inputPSD;
  end
  
  currentPSDLength  = numel(currentFreqAxes);
  currentXAxis      = repmat(xAxisPosition(idx), [1, currentPSDLength]);

  plot3(hAx, currentXAxis, currentFreqAxes, currentFreqValues, 'r', 'LineWidth', 1.0);
  hold(hAx, 'on')

end

% 
% currentXTick       = get(hAx, 'XTick');
% currentXTickLabels = get(hAx, 'XTickLabels');
% 
% set(hAx, 'XTickLabels', xAxisTickName);