% define input signals

% lisftOfSignals = {'M1_MINUS_M2', 'M3_MINUS_M4', 'R1_MINUS_R2'};
lisftOfSignals = {'M1_MINUS_M2', 'M3_MINUS_M4'};
variables.fs   = 1e3;

path2CartoFile     = fullfile('data','Auriculas','#4','Carto1','EG.mat');
fullPath2CartoFile = fullfile(rootFolder(), path2CartoFile); 

listOfPoints = 'all';
[data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile);

allPoints = unique(pointsInfo);  
numPoints  = numel(allPoints);

useSignalFlag = char(1, numPoints);
%% loop through the adquisition points

for idxPoint = 1:10
%for idxPoint=1:numPoints
    
currentResults = [];
currentPoint   = allPoints(idxPoint);
  
idx      = find(pointsInfo == currentPoint);
ecgData  = data(idx);
    

fprintf('\nProcessing point %d\n', currentPoint);
   

%% Check all signals are present and preprocess input data.
dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;
dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;

[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);
assert(all(signalPresent), 'Missing signals');

selDataIn.signalData  = dataIn.signalData(idxSignals, :);
selDataIn.signalNames = dataIn.signalNames(idxSignals)  ;
selDataIn.fs          = dataIn.fs;

% Preprocess
% Schilling (p145) sets filters to (30, 150)
[numSignals, lengthSignal] = size(selDataIn.signalData);

subPlotIdx = 0;
figure();

for currentIdx = 1:numSignals
    
    currentSignal     = selDataIn.signalData(currentIdx, :);
    currentSignalName = selDataIn.signalNames{currentIdx};
    
      % get the PSD for the current input signal
    [inputPSD, inputPSDOmega] = ...
      getPowerSpectralDensity(currentSignal, selDataIn.fs, ...
      'windowType', 'rectangular'        , ...
      'windowSize', length(currentSignal), ...
      'lengthFFT' , 8192, 'numberOverlap', [], 'psdCutOff', 180);
  
    % Time domain 
    subPlotIdx = 1 + subPlotIdx;
    subplot(numSignals, 2, subPlotIdx),
    plot(currentSignal, 'Color', [1 0 0], 'LineWidth', 1.5);
    signalTitle = sprintf('%s - point %d', currentSignalName, currentPoint);
    title(signalTitle, 'Interpreter', 'none');
  
    
    % Freq domain
    subPlotIdx = 1 + subPlotIdx;
    subplot(numSignals, 2, subPlotIdx),
    plot(inputPSDOmega, inputPSD, 'LineWidth', 2);
    
%     % get the axes
% %     f = gcf();
% %     figAxes = gca();
% 
%     title(f.Children(2), sprintf('Point %d', currentSignalName, currentPoint));
%     
%     subplot(figAxes, 1,2,1)
%     
%     figHandle = figure();
%     subplot(figHandle.CurrentAxes, 1,2,1)
%     
    
end
%     prompt = 'Keep this signal? Y/N [Y]: ';
%     useSignalFlag(idxPoint) = input(prompt,'s');
    
    pause(1);
    close(gcf);

end

%         
% f1HighPassCutOff = 1.0;
% f2StopBandCutOff = 50 ;
% f3LowPassCutOff  = 250;
% dataForICA = preprocessForICAv2(selDataIn, ...
%     f1HighPassCutOff, f2StopBandCutOff, f3LowPassCutOff);
% 
% 
% plotSignalAndSpectrum(dataForICA(idx, :), samplingFreq);