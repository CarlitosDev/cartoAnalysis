function plotDataAnalysisParallelAndMeridianV2(cartoPointAnalysis, meshData, xData, yData , zData, selectedPointsData)

%% PLOTDATAANALYSISPARALLELANDMERIDIANV2 plot the results of the parallel/meridian analysis
%
%
% Carlos Aguilar  09/11/2016
% TO-DO: Sort out the pointsIdx and pointsId malarkey

%% hard-coded parameters

    numSamples        = selectedPointsData.lengthSignal;
    samplingFrequency = selectedPointsData.fs;
    timeLength        = numSamples/samplingFrequency;

    doScaling         = true;
    xAxisStretch      = 3;
    
    doPlotMesh        = true;
    
% define colours    
    meshColour         = [232,255,255]/255;
    markersColour      = [1 0 0];
    parallelLineColour = [255,127,80] /255;
    meridianLineColour = [255 140 0]  /255;
    allPointsColour    = [240, 240, 240]/255;

%% inliners
  
    scaleFrom0To1      = @(x) (x-min(x))/(max(x)-min(x));
    vectorToStr        = @(x) cellstr(num2str(x));
    
%% Input data

  selectedPointIdx = selectedPointsData.selectedPointIdx;
  
  parallelPath = selectedPointsData.parallelPath;
  meridianPath = selectedPointsData.meridianPath;
  pathInfo     = selectedPointsData.pathInfo;
  
  % bundle all the points to plot them as small spheres and add the point
  % number
  idxAllPoints      = unique([selectedPointsData.allParallelPointsIdx , ...
                              selectedPointsData.allMeridianPointsIdx]);
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
    text(selectedPointsData.textPosition(idxAllPoints(i), 1), ...
         selectedPointsData.textPosition(idxAllPoints(i), 2), ...
         selectedPointsData.textPosition(idxAllPoints(i), 3), ...
         sprintf('%d', selectedPointsData.pointsId(idxAllPoints(i))));
  end
  
  
  set(0, 'DefaultTextInterpreter', 'none');
  
  xLabelString = ...
    sprintf('Signal %s\n%s\n Point selected %d\n', ...
    selectedPointsData.inputSignalName, ...
    selectedPointsData.testType, ...
    selectedPointsData.selectedPoint);  
  xlabel(xLabelString);
 
%% plot the 3D mesh, the parallels and meridians and the selected points

  % plot all the points belonging to the parallels and meridians
  plot3(xData,yData, zData, ...
        'o', 'MarkerSize' , 4, ...
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
        'LineWidth', 1.5 );

% plot parallel points      
  plot3(xData(pathInfo.allParallelPointsIdx) , ...
        yData(pathInfo.allParallelPointsIdx) , ...
        zData(pathInfo.allParallelPointsIdx) , ...
        'MarkerFaceColor', parallelLineColour, ...
        'LineWidth', 1.5 );
  
% point the plot to the selected point
  viewAz = rad2deg(selectedPointsData.selectedAz);
  viewEl = rad2deg(selectedPointsData.selectedEl);
  cAxes  = gca();
  cAxes.View = [viewAz, viewEl];
  axis(cAxes, 'tight')
  
  drawnow
  
%% Get the point Id's

    % parallel
    parallelPointsList    = ...
      selectedPointsData.pointsId(pathInfo.allParallelPointsIdx);
    xAxisParallelTickName = vectorToStr(parallelPointsList');
    
    % meridian
    meridianPointsList    = ... 
      selectedPointsData.pointsId(pathInfo.allMeridianPointsIdx);
    xAxisMeridianTickName = vectorToStr(meridianPointsList');
    
  
%% (2.1) Time-Plot EGM parallel
  
  hAxTime  = subplot(2,3,[2,2]);
  timeAxes = linspace(0, timeLength, numSamples);
  cumDistParallel = [0, cumsum([parallelPath.distAB])];

  % first plot the selected point
    idx = 1;
    idxCurrentParallel = parallelPath(idx).idxPointA;

    % build a x-axis to place the signal
    currentXAxis       = zeros([1, numSamples]);
    currentCartoData   = cartoPointAnalysis{idxCurrentParallel};
    if doScaling
      currentTimeValues = scaleFrom0To1(currentCartoData.dataForICA);
    else
      currentTimeValues = currentCartoData.dataForICA; %#ok<UNRCH>
    end

    plot3(hAxTime, currentXAxis, ...
        timeAxes, currentTimeValues, 'r', ...
        'LineWidth', 1.2);
    hold(hAxTime, 'on');      
    
  % plot the remaining signals
  for idx = 1:pathInfo.numParallelPoints-1

    idxCurrentParallel = parallelPath(idx).idxPointB;
    
    % shift the axis by the distance between points
    currentXAxis = currentXAxis + parallelPath(idx).distAB;
  
    currentCartoData = cartoPointAnalysis{idxCurrentParallel};
    if doScaling
        currentTimeValues = scaleFrom0To1(currentCartoData.dataForICA);
    else
        currentTimeValues = currentCartoData.dataForICA; %#ok<UNRCH>
    end

    plot3(hAxTime , currentXAxis          , ...
          timeAxes, currentTimeValues, 'b', ...
          'LineWidth', 1.0);
   
    hold(hAxTime, 'on');

  end
  
  hAxTime.View = [-65, 80];
  axis(hAxTime, 'tight');
  axis(hAxTime, 'ij');
  
  % set point-numbers
  hAxTime.XTick      = cumDistParallel;
  hAxTime.XTickLabel = xAxisParallelTickName;
  hAxTime.FontAngle  = 'italic';
  hAxTime.FontSize   = 7;
  grid(hAxTime, 'on');
  parallelTitle = sprintf('%d Points along the parallel', ...
    selectedPointsData.numParallelPoints);
  title(hAxTime , parallelTitle)
  ylabel(hAxTime, 't(seconds)');


%% (2.2)  Plot frequency EGM parallel

  hAx = subplot(2,3,[5,5]);

  % first plot the selected point
    idx = 1;
    idxCurrentParallel = parallelPath(idx).idxPointA;

    currentCartoData   = cartoPointAnalysis{idxCurrentParallel};
    currentFreqAxes    = currentCartoData.inputPSDOmega;

    currentPSDLength   = numel(currentFreqAxes);
    currentXAxis       = zeros([1, currentPSDLength]);

    if doScaling
        currentFreqValues = scaleFrom0To1(currentCartoData.inputPSD);
    else
        currentFreqValues = currentCartoData.inputPSD; %#ok<UNRCH>
    end
    plot3(hAx, currentXAxis, ...
        currentFreqAxes, currentFreqValues, ...
        'r', 'LineWidth', 1.0);
    hold(hAx, 'on');
  
  for idx = 1:pathInfo.numParallelPoints-1

       idxCurrentParallel = parallelPath(idx).idxPointB;
      
      currentCartoData  = cartoPointAnalysis{idxCurrentParallel};
      currentFreqAxes   = currentCartoData.inputPSDOmega;
      if doScaling
          currentFreqValues = scaleFrom0To1(currentCartoData.inputPSD);
      else
          currentFreqValues = currentCartoData.inputPSD; %#ok<UNRCH>
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
  
    hAxMrdn  = subplot(2,3,[3,3]);
    timeAxes = linspace(0, timeLength, numSamples);
    cumDistMeridian = [0, cumsum([meridianPath.distAB])];

    % first plot the selected point
    idx = 1;
    idxCurrentMeridian = meridianPath(idx).idxPointA;

    % build a x-axis to place the signal
    currentXAxis       = zeros([1, numSamples]);
    currentCartoData   = cartoPointAnalysis{idxCurrentMeridian};
    if doScaling
        currentTimeValues = scaleFrom0To1(currentCartoData.dataForICA);
    else
        currentTimeValues = currentCartoData.dataForICA; %#ok<UNRCH>
    end

    plot3(hAxMrdn, currentXAxis, ...
        timeAxes, currentTimeValues, 'r', ...
        'LineWidth', 1.2);
    
    hold(hAxMrdn, 'on');      

    for idx = 1:pathInfo.numMeridianPoints-1

        idxCurrentMeridian = meridianPath(idx).idxPointB;

        % shift the axis by the distance between points
        currentXAxis = currentXAxis + meridianPath(idx).distAB;

        currentCartoData = cartoPointAnalysis{idxCurrentMeridian};
        if doScaling
            currentTimeValues = scaleFrom0To1(currentCartoData.dataForICA);
        else
            currentTimeValues = currentCartoData.dataForICA; %#ok<UNRCH>
        end

        plot3(hAxMrdn , currentXAxis          , ...
            timeAxes, currentTimeValues, 'b', ...
            'LineWidth', 1.0);

        hold(hAxMrdn, 'on');

    end

    hAxMrdn.View = [-65, 80];
    axis(hAxMrdn, 'tight');
    axis(hAxMrdn, 'ij');

    % set point-numbers
    hAxMrdn.XTick      = cumDistMeridian;
    hAxMrdn.XTickLabel = xAxisMeridianTickName;
    hAxMrdn.FontAngle  = 'italic';
    hAxMrdn.FontSize   = 7;
    grid(hAxMrdn, 'on');
    meridianTitle = sprintf('%d Points along the meridian', ...
    selectedPointsData.numMeridianPoints);
    title(hAxMrdn , meridianTitle)    
    ylabel(hAxMrdn, 't(seconds)');
  
  
%% (3.2) Plot frequency EGM parallel
  hAxMrdnFreq = subplot(2,3,[6,6]);

  % first plot the selected point
  idx = 1;
  idxCurrentMeridian = meridianPath(idx).idxPointA;
  
  currentCartoData   = cartoPointAnalysis{idxCurrentMeridian};
  currentFreqAxes    = currentCartoData.inputPSDOmega;
  
  currentPSDLength   = numel(currentFreqAxes);
  currentXAxis       = zeros([1, currentPSDLength]);
  
  if doScaling
      currentFreqValues = scaleFrom0To1(currentCartoData.inputPSD);
  else
      currentFreqValues = currentCartoData.inputPSD; %#ok<UNRCH>
  end
  plot3(hAxMrdnFreq, currentXAxis, ...
      currentFreqAxes, currentFreqValues, ...
      'r', 'LineWidth', 1.0);
    
	hold(hAxMrdnFreq, 'on');
  
  for idx = 1:pathInfo.numMeridianPoints-1

      idxCurrentMeridian = meridianPath(idx).idxPointB;
      
      currentCartoData  = cartoPointAnalysis{idxCurrentMeridian};
      currentFreqAxes   = currentCartoData.inputPSDOmega;
      if doScaling
          currentFreqValues = scaleFrom0To1(currentCartoData.inputPSD);
      else
          currentFreqValues = currentCartoData.inputPSD; %#ok<UNRCH>
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
  
  drawnow