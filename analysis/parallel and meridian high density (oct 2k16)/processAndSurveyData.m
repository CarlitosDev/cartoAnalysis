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

  
  % Preprocessing features
  f1HighPassCutOff = 1.0;
  f2StopBandCutOff = 50 ; % notch
  f3LowPassCutOff  = 150;
  
  
%% Pick a point to research.
  
currentPoint = 886;
idx = pointsInfo.pointToIdxMap(currentPoint);
ecgData  = data(idx);

dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;


signalData  = cell2mat(ecgData.signal');


% set the signals to analyse
pentaArray = listOfPentaArraySignals();

signalToPlot = pentaArray.pentaBipolar{4};
[~, refSignalIdx] = ismember(signalToPlot, ecgData.tipo_ECG);
referenceSignal = signalData(refSignalIdx, :);
  set(0, 'DefaultTextInterpreter', 'none');
  
  figTitle = sprintf('Point%d_%s', currentPoint, signalToPlot);
  
  plotSignalAndSpectrum(referenceSignal, 1e3);
  title(figTitle);
  
  
  

 saveFigureAsPNG(figTitle, rootFolder);
% 
%  close(gcf());
  
  
  %%
  % Schilling (p145) sets filters to (30, 150)
  dataIn = [];
  dataIn.signalData = signalData;
  dataIn.fs = 1e3;
  
  [numSignals, lengthSignal] = size(signalData);
  
  dataProcessed = preprocessForICAv2(dataIn, ...
    f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);
  
  % sample the original signal for the b-spline in chunks of 20% the total size
  samplingStep = 20;
  fprintf('Removing baseline wandering...');
  for idx=1:numSignals
    [dataProcessed(idx, :), ~] = ...
      removeBaselineWanderingV3(dataProcessed(idx, :), samplingStep);
  end
  fprintf('done\n');
  
  %%
  
  
  %%
%   [powSpectDens, psdOmega] = getPowerSpectralDensity(signalData(1, :), dataIn.fs);  
% figure, plot(psdOmega, powSpectDens);
  %%
  
% idxSignal = 20;
%   
% figure,
% currentSignal = signalData(idxSignal, :);
% plot(currentSignal), hold on;
% title('rawdata');
% currentSignal = dataProcessed(idxSignal, :);
% plot(currentSignal), hold on;
% title('dataProcessed');


  
  
  %%
  dataIn = [];
  dataIn.signalNames = ecgData.tipo_ECG;
  dataIn.signalData = dataProcessed;

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
% fastICA = [];
% fastICA.approach = 'defl'; %'symm'; % 'defl'
% fastICA.g = 'gauss'; %'tanh', 'gauss', 'skew','pow3'
% fastICA.stabilization = 'on';
% numberOfICs = 2;
% 
%   [signalICA, A, W] = fastica(dataIn.signalData, 'approach', fastICA.approach, ...
%        'g', fastICA.g, 'stabilization', fastICA.stabilization, 'numOfIC', numberOfICs);
% 
%      figure,
%     for idxSignal=1:numberOfICs
%       
%       idxInternalSignal = idxSignal+(idx-1)*signalsPerBatch;
%       plot(dataIn.signalData(idxInternalSignal, :));
%       hold on
%       title(sprintf('IC %d', idxSignal));
%       pause
%       
%     end     
