function plotDataAnalysisParallelAndMeridian(cartoPointAnalysis, meshData, selectedPointsData)


%% hard-coded parameters (review this...)

    numSamples        = 2500;
    samplingFrequency = 1e3;
    timeLength        = numSamples/samplingFrequency;

    doScaling         = true;
    xAxisStretch      = 3;

%% Input data

  selectedPointIdx  = selectedPointsData.selectedIdx;
  idxAllPoints      = selectedPointsData.idxAllPoints;
  numSelectedPoints = numel(idxAllPoints);
  
  idxParallel       = selectedPointsData.idxParallel;
  idxMeridian       = selectedPointsData.idxMeridian;


  %% inliners

    scaleFrom0To1      = @(x) (x-min(x))/(max(x)-min(x));
    vectorToStr        = @(x) cellstr(num2str(x));

%% Plot mesh and selected points

  h = figure;
  set(h, 'units', 'normalized', 'outerposition', [0 0 1 1]);

  subplot(2,3,[4,1]);
  
	% plot mesh
  trimesh(meshData.faces, ...
          meshData.X, meshData.Y, meshData.Z, ...
          'FaceColor', 'none', ...
          'EdgeColor', 220*ones(1,3)/255);
  
  % Add points Id  
  for i=1:numSelectedPoints
    hold on;
    text(selectedPointsData.textPosition(i, 1), ...
         selectedPointsData.textPosition(i, 2), ...
         selectedPointsData.textPosition(i, 3), ...
         sprintf('%d', selectedPointsData.pointsId(i)));
  end
    
  xLabelString = ...
    sprintf('Point selected %d\nParallel %d points. Meridian %d points', ...
    selectedPointIdx, numel(idxParallel), numel(idxMeridian));
  xlabel(xLabelString);
 
 
  % plot parallel points
  plot3(selectedPointsData.xDataParallel, ...
        selectedPointsData.yDataParallel, ...
        selectedPointsData.zDataParallel, ...
        'o', 'MarkerSize' , 10, ...
        'MarkerFaceColor',[0 1 1]);
      
  plot3(selectedPointsData.xDataParallel, ...
        selectedPointsData.yDataParallel, ...
        selectedPointsData.zDataParallel, ...
        'MarkerFaceColor',[255,127,80]/255, ...
        'LineWidth', 1.5 );

  % plot meridian points
  plot3(selectedPointsData.xDataMeridian, ...
        selectedPointsData.yDataMeridian, ...
        selectedPointsData.zDataMeridian, ...
        'o', 'MarkerSize' , 10, ...
        'MarkerFaceColor',[255 140 0]/255);   
      
  plot3(selectedPointsData.xDataMeridian, ...
        selectedPointsData.yDataMeridian, ...
        selectedPointsData.zDataMeridian, ...
        'MarkerFaceColor',[1 0 0], ...
        'LineWidth', 1.5 );
      
  % plot selected point
  plot3(selectedPointsData.xDataMeridian(1), ...
        selectedPointsData.yDataMeridian(1), ...
        selectedPointsData.zDataMeridian(1), ...
        'o', 'MarkerSize' , 20, ...
        'MarkerFaceColor',[1 1 0]);
      
  viewAz = rad2deg(selectedPointsData.selectedAz);
  viewEl = rad2deg(selectedPointsData.selectedEl);
  cAxes  = gca();
  cAxes.View = [viewAz, viewEl];
  axis(cAxes, 'tight')

%% Distances along paths
% Get distances between points of a path. Use Euclidean for the time being

  % PARALLEL distances
  allParallel = horzcat(selectedPointsData.xDataParallel, ...
                        selectedPointsData.yDataParallel, ...
                        selectedPointsData.zDataParallel);

  % distances from pointA to point B. Take the distance matrix and arrange it
  pDist = squareform(pdist(allParallel, 'euclidean'));
  [numRows, numCols] = size(pDist);
  linearInd     = sub2ind([numRows, numCols], 1:numRows-1, 2:numCols);
  parallelDist  = pDist(linearInd);
  xAxisParallel = [0, cumsum(parallelDist)];
  

  % MERIDIAN distances
  allMeridian = horzcat(selectedPointsData.xDataMeridian, ...
                        selectedPointsData.yDataMeridian, ...
                        selectedPointsData.zDataMeridian);

  % distances from pointA to point B. Take the matrix and a
  pDist = squareform(pdist(allMeridian, 'euclidean'));
  [numRows, numCols] = size(pDist);
  linearInd     = sub2ind([numRows, numCols], 1:numRows-1, 2:numCols);
  meridianDist  = pDist(linearInd);
  xAxisMeridian = [0, cumsum(meridianDist)];
  
  %% Get the point Id's
  
  % parallel
  [~,idxInA,~]   = intersect(idxAllPoints, idxParallel);
  parallelPoints = selectedPointsData.pointsId(idxInA)';
  xAxisParallelTickName = vectorToStr(parallelPoints);
  
  % meridian
  [~,idxInA,~]   = intersect(idxAllPoints, idxMeridian);
  meridianPoints = selectedPointsData.pointsId(idxInA)';
  xAxisMeridianTickName = vectorToStr(meridianPoints);


%% (2.1) Time-Plot EGM parallel

  hAxTime  = subplot(2,3,[2,2]);
  timeAxes = linspace(0, timeLength, numSamples);
  currentXAxisBase = ones([1, numSamples])*xAxisStretch;

  for idx = 1:numel(idxParallel)

    currentCartoData  = cartoPointAnalysis{idxParallel(idx)};
    if doScaling
        currentTimeValues = scaleFrom0To1(currentCartoData.dataForICA);
    else
        currentTimeValues = currentCartoData.dataForICA; %#ok<UNRCH>
    end

    currentXAxis = xAxisParallel(idx)*currentXAxisBase;

    if idxParallel(idx) == selectedPointIdx
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
  
  hAxTime.View = [-65, 80];
  axis(hAxTime, 'tight');
  axis(hAxTime, 'ij');
  
  % set point-numbers
  hAxTime.XTick      = xAxisParallel*xAxisStretch;
  % using the distances is very confusing
  %hAxTime.XTickLabel = cellstr(num2str(xAxisParallel', '%3.2f'));
  hAxTime.XTickLabel = xAxisParallelTickName;
  hAxTime.FontAngle  = 'italic';
  hAxTime.FontSize   = 7;
  grid(hAxTime, 'on');
  title(hAxTime , 'Points along the parallel')
  ylabel(hAxTime, 't(seconds)');
  

%% (2.2)  Plot frequency EGM parallel

  hAx = subplot(2,3,[5,5]);

  for idx = 1:numel(idxParallel)

    currentCartoData  = cartoPointAnalysis{idxParallel(idx)};
    currentFreqAxes   = currentCartoData.inputPSDOmega;
    if doScaling
        currentFreqValues = scaleFrom0To1(currentCartoData.inputPSD);
    else
        currentFreqValues = currentCartoData.inputPSD; %#ok<UNRCH>
    end

    currentPSDLength = numel(currentFreqAxes);
    currentXAxis     = xAxisParallel(idx)*ones([1, currentPSDLength])*xAxisStretch;

      if idxParallel(idx) == selectedPointIdx
        plot3(hAx, currentXAxis, ...
        currentFreqAxes, currentFreqValues, ...
        'r', 'LineWidth', 1.0);
      else
        plot3(hAx, currentXAxis, ...
        currentFreqAxes, currentFreqValues, ...
        'b', 'LineWidth', 1.0);
      end

    hold(hAx, 'on');

  end

  axis(hAx, 'tight')
  axis(hAx, 'ij');
  hAx.View = [-65, 80];
  hAx.XTick      = xAxisParallel*xAxisStretch;
  hAx.XTickLabel = xAxisParallelTickName; 
  hAx.FontAngle  = 'italic';
  hAx.FontSize   = 8;
  grid(hAx, 'on');
  ylabel(hAx, 'f(Hertz)');
  
  
  %% (3.1)  Time-Plot EGM meridian

  hAxMrdn  = subplot(2,3,[3,3]);
  timeAxes = linspace(0, timeLength, numSamples);
  currentXAxisBase = ones([1, numSamples])*xAxisStretch;

  for idx = 1:numel(idxMeridian)

    currentCartoData  = cartoPointAnalysis{idxMeridian(idx)};
    if doScaling
        currentTimeValues = scaleFrom0To1(currentCartoData.dataForICA);
    else
        currentTimeValues = currentCartoData.dataForICA; %#ok<UNRCH>
    end

    currentXAxis = xAxisMeridian(idx)*currentXAxisBase;

    if idxMeridian(idx) == selectedPointIdx
    plot3(hAxMrdn, currentXAxis, ...
          timeAxes, currentTimeValues, 'r', ...
          'LineWidth', 1.2);
    else
    plot3(hAxMrdn, currentXAxis, ...
          timeAxes, currentTimeValues, 'b', ...
          'LineWidth', 1.0);
    end
    hold(hAxMrdn, 'on');

  end
  
  hAxMrdn.View = [-65, 80];
  axis(hAxMrdn, 'tight')
  axis(hAxMrdn, 'ij');
  
  % set point-numbers
  hAxMrdn.XTick      = xAxisMeridian*xAxisStretch;
  % using the distances is very confusing
  %hAxMrdn.XTickLabel = cellstr(num2str(xAxisMeridian', '%3.2f'));
  hAxMrdn.XTickLabel = xAxisMeridianTickName;
  hAxMrdn.FontAngle  = 'italic';
  hAxMrdn.FontSize   = 7;
  grid(hAxMrdn, 'on');
  title(hAxMrdn, 'Points along the meridian');
  ylabel(hAxMrdn, 't(seconds)');
  
  
%% (3.2) Plot frequency EGM parallel
hAxMrdnFreq = subplot(2,3,[6,6]);

  for idx = 1:numel(idxMeridian)

    currentCartoData  = cartoPointAnalysis{idxMeridian(idx)};
    currentFreqAxes   = currentCartoData.inputPSDOmega;
    if doScaling
        currentFreqValues = scaleFrom0To1(currentCartoData.inputPSD);
    else
        currentFreqValues = currentCartoData.inputPSD; %#ok<UNRCH>
    end

    currentPSDLength = numel(currentFreqAxes);
    currentXAxis     = xAxisMeridian(idx)*ones([1, currentPSDLength])*xAxisStretch;  

      if idxMeridian(idx) == selectedPointIdx
        plot3(hAxMrdnFreq, currentXAxis, ...
        currentFreqAxes, currentFreqValues, ...
        'r', 'LineWidth', 1.0);
      else
        plot3(hAxMrdnFreq, currentXAxis, ...
        currentFreqAxes, currentFreqValues, ...
        'b', 'LineWidth', 1.0);
      end

    hold(hAxMrdnFreq, 'on');

  end

  axis(hAxMrdnFreq, 'tight')
  axis(hAxMrdnFreq, 'ij');
  hAxMrdnFreq.View = [-65, 80];
  hAxMrdnFreq.XTick      = xAxisMeridian*xAxisStretch;
  hAxMrdnFreq.XTickLabel = xAxisMeridianTickName; 
  hAxMrdnFreq.FontAngle  = 'italic';
  hAxMrdnFreq.FontSize   = 8;
  grid(hAxMrdnFreq, 'on');
  ylabel(hAxMrdnFreq, 'f(Hertz)');