function plotCatheterSignalsFrequencyV2(cathTable, currentPoint, currentType)

% PLOTCATHETERSIGNALSFREQUENCYV2 plot Pentaray spatial, temporal and
% frequencial colouring by branch of the catheter.
%
% cathTable: A table containing the catheter points. Can be generated with
% 'loadRamonYCajalData.m'
% 
% currentPoint: A valid pointId from the Pentaray list
% 
% currentType: Either 'MON' or 'BIP'
%
% NOTE: same as plotCatheterSignalsFrequency but colouring the branches.



%%

  ancho = 1.8; % for the OI
   
  idxSubSet  = cathTable{:, 'pentaRayPoint'} == currentPoint;
  subSet     = cathTable(idxSubSet, :);
  idxMon     = strcmpi(subSet{:, 'type'}, currentType);
  subSet     = subSet(idxMon, :);
  numSignals = height(subSet);  
  
  cathCenter = mean([subSet{:, {'x', 'y', 'z'}}]);

  hFig = figureMax();
  lAx  = subplot(2,2,1);
  hold(lAx, 'on');
  
	tableAx = subplot(2,2,3);

  plot3(lAx, subSet{:, 'x'}, ...
    subSet{:, 'y'}, subSet{:, 'z'}, 'Marker', 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', [255,127,80]./255, ...
    'LineStyle', 'none');
  
	plot3(lAx, cathCenter(1), cathCenter(2), cathCenter(3), ...
  'MarkerSize', 15, 'MarkerFaceColor', [1 0 0.2], 'Marker', 'o');

  lAx.View = [-10, 30];
  title(lAx, ...
    sprintf('%d - %s catheter placement', currentPoint, currentType));

  cMap = colormap('lines');
  
  % Add the catherer lines
  if strcmpi(currentType, 'BIP')
    intA = 1:3:numSignals;
    intB = [intA(2:end)-1, numSignals];
    fCutOff = 100;
    numElec = 3;
  else
    intA = 1:4:numSignals;
    intB = [intA(2:end)-1, numSignals];
    fCutOff = 30;
    numElec = 4;
  end
  
  for iChunk = 1:numel(intA)
    x = [subSet{intA(iChunk):intB(iChunk), 'x'}; cathCenter(1)];
    y = [subSet{intA(iChunk):intB(iChunk), 'y'}; cathCenter(2)];
    z = [subSet{intA(iChunk):intB(iChunk), 'z'}; cathCenter(3)];   
    plot3(lAx,x,y,z, 'Color', cMap(iChunk, :), 'LineWidth', 2);
  end
  
  
  rAx1 = subplot(2,2,2);
  rAx2 = subplot(2,2,4);
  hold(rAx1, 'on');
  hold(rAx2, 'on');

  currentSignal = subSet{1, 'signal'};
  numSamples    = length(currentSignal);
  samplingFreq  = 1e3;
  timeLength    = numSamples/samplingFreq;
  timeAxes      = linspace(0, timeLength, numSamples);

  rAx1.View = [-75.26, 62.64];
  rAx2.View = [104,38];
 
  parametersTemp = [];
  idxColour      = 0;
  
  for idx=1:numSignals
    
    if mod(idx-1, numElec) == 0
      idxColour = idxColour + 1;
    end
    
    currentSignal = subSet{idx, 'signal'};
    
    text(lAx, 0.3+subSet{idx, 'x'}, ...
      subSet{idx, 'y'}, subSet{idx, 'z'}, ...
      subSet{idx, 'pentaRayElectrode'}, 'Interpreter', 'none');
    
    % signal - timeplot
    currentXAxis = repmat(idx, [1, numSamples]);
    
    plot3(rAx1, currentXAxis, timeAxes, currentSignal, ...
      'b', 'LineWidth', 1.2, 'Color', cMap(idxColour, :));
    
    % signal - freqplot
    % get the PSD for the input signal
    [currentFreqValues, currentFreqAxes] = ...
      getPowerSpectralDensity(currentSignal, samplingFreq, ...
      'windowType', 'rectangular', 'windowSize', numSamples, ...
      'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', fCutOff);
    [~, idxPeak] = max(currentFreqValues);
    maxFrequency = currentFreqAxes(idxPeak);
    currentPSDLength  = numel(currentFreqAxes);
    currentXAxis      = repmat(idx, [1, currentPSDLength]);
    
    plot3(rAx2, currentXAxis, currentFreqAxes, ...
      currentFreqValues, 'r', 'LineWidth', 1.0, ...
      'Color', cMap(idxColour, :));

  % get some parameters
    [botFreq,~,~,~,periodogramInput,fz] = df_Ng(currentSignal, samplingFreq);
    oi = organizationIndex(periodogramInput,fz,botFreq,ancho);
    
    parametersTemp(idx).signalId = subSet{idx, 'pentaRayElectrode'};
    parametersTemp(idx).botFreq  = botFreq;
    parametersTemp(idx).fPeak    = maxFrequency;
    parametersTemp(idx).OI       = oi;
    parametersTemp(idx).meanVol  = mean(currentSignal);
    parametersTemp(idx).stdVol   = std(currentSignal);
    parametersTemp(idx).eKurt    = kurtosis(currentSignal)-3.0;

  end
  
  ylabel(rAx1, 'Time(s)')
  rAx1.View = [-75.26, 62.64];
	ylabel(rAx2, 'Freq(Hz)')
  rAx2.View = [104,38];
  rAx1.XDir = 'normal';
  rAx1.YDir = 'reverse';
  rAx2.YDir = 'normal';
	rAx2.XDir = 'reverse';
  
  rAx1.XTick      = 1:numSignals;
  rAx1.XTickLabel = subSet{:, 'pentaRayElectrode'};
  set(rAx1, 'TickLabelInterpreter', 'none');
  
	rAx2.XTick      = 1:numSignals;
  rAx2.XTickLabel = subSet{:, 'pentaRayElectrode'};
  set(rAx2, 'TickLabelInterpreter', 'none');
  
	parametersTable = struct2table(parametersTemp);
  delete(tableAx);

  listOfVars = {'botFreq', 'fPeak', 'OI', 'meanVol', 'stdVol', 'eKurt'};
  uiAxes = uitable(hFig, 'Data', parametersTable{:, listOfVars}, ...
    'ColumnName',listOfVars, ...
    'RowName', parametersTable{:, 'signalId'});
  
  rAx2.Units      = 'pixels';
  currentPos      = uiAxes.Position;
  currentPos(3)   = rAx2.Position(1) - 100;
  uiAxes.Position = currentPos;
  
  

