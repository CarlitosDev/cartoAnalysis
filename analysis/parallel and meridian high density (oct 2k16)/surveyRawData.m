%% survey data
%   Load mat file and dump into workspace

  rawDataFileName = 'rawData_FA_RamonYCajal_ECG.mat';

  matchingFNames = which(rawDataFileName);
  if isempty(matchingFNames)
    errordlg('Can''t resolve path to file. Run initialise()');
    return;
  end

  rawData = load(matchingFNames);

  fNames = fieldnames(rawData);
  for idx=1:numel(fNames) 
    assignin('base', fNames{idx},  rawData.(fNames{idx}));
  end

%% Pick a point to research.
  
  currentPoint = pointsInfo(100);
  
  idx      = find(pointsInfo == currentPoint);
  ecgData  = data(idx);
  
  rawData  = cell2mat(ecgData.signal');
  
  %basicStats = Utils_basic_stats( rawData );
  
  [inputDataSorted, meanEnergyValues, idxNRG] = sortDataByMeanEnergy(rawData);

  dataIn = [];
  dataIn.signalNames = ecgData.tipo_ECG(idxNRG);
  dataIn.signalData  = inputDataSorted;
  
  numSignals = numel(dataIn.signalNames);

%% Box plot of the signals 
% Split into 4 frames
  
  figureLeft();
  numBatches = 4;
  signalsPerBatch = numSignals/numBatches;
    
  for idx=1:numBatches
    subplot(2,2,idx),
    sliceBounds = [1+(idx-1)*signalsPerBatch, idx*signalsPerBatch];
    boxplot(dataIn.signalData(sliceBounds(1):sliceBounds(2), :)', ...
            'orientation','horizontal', ...
            'labels', dataIn.signalNames(sliceBounds(1):sliceBounds(2)))
    title(sprintf('Boxplot for point %d (%d/%d)', currentPoint, idx, numBatches));          
  end

%% Time plot of the signals 
% Split into 4 frames and add some offset to ease the viz.
  
  prevOffset = 0;
  figure
  for idx=1:numBatches
    subplot(2,2,idx)
    for idxSignal=1:signalsPerBatch
      
      idxInternalSignal = idxSignal+(idx-1)*signalsPerBatch;
      currentSignal     = prevOffset + ...
        dataIn.signalData(idxInternalSignal, :);
      plot(currentSignal);
      hold on
      title(dataIn.signalNames{idxInternalSignal});
      pause
      
      prevOffset = prevOffset + ...
        max(abs(dataIn.signalData(idxInternalSignal, :)));
    end
  end
  
 %% Set fastICA parameters in here
fastICA = [];
fastICA.approach = 'defl'; %'symm'; % 'defl'
fastICA.g = 'gauss'; %'tanh', 'gauss', 'skew','pow3'
fastICA.stabilization = 'on';
numberOfICs = 2;

  [signalICA, A, W] = fastica(dataIn.signalData, 'approach', fastICA.approach, ...
       'g', fastICA.g, 'stabilization', fastICA.stabilization, 'numOfIC', numberOfICs);

     figure,
    for idxSignal=1:numberOfICs
      
      idxInternalSignal = idxSignal+(idx-1)*signalsPerBatch;
      plot(dataIn.signalData(idxInternalSignal, :));
      hold on
      title(sprintf('IC %d', idxSignal));
      pause
      
    end     
