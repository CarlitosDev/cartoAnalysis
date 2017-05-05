%%
%    Script to generate synth data
    
rootFolder  = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos';
addpath(genpath(rootFolder));
mainFolder  = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos\tests\filter benchmark\synthesized data';
currentTest = 'test1';


%% Synth ECG + noise
% Call ecgsyn
% IEEE Transactions On Biomedical Engineering, 50(3), 289-294, March 2003.
samplingFreq = 1e3;
numBeats = 12;
Anoise = 0;
hrmean = 72;
hrstd = 1;
lfhfratio = 0.5;
sfint = 1000;
ti=[-70 -15 0 15 100];
ai=[1.2 -5 30 -7.5 0.75];
bi=[0.25 0.1 0.1 0.1 0.4];

[simECG, ipeaks, addNoise] = ecgsyn(samplingFreq,numBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);
cleanECG = simECG - addNoise;

% not very glad on using previous noise...
simECG = cleanECG;

%whiteNoise = awgn(cleanECG, 10, 'measured'); % Add white Gaussian noise.


%% Synth AF

numSamples = length(simECG);
atrialFibrillation = synthesizeAtrialFibrillation(numSamples, samplingFreq)  ;

% scale AF(uV) to mV
atrialFibrillation = atrialFibrillation/1e3;


%% check out input values. Scale signals...

statsCleanECG = Utils_basic_stats(cleanECG);
statsNoiseECG = Utils_basic_stats(addNoise);
statsAF       = Utils_basic_stats(atrialFibrillation);

scale0To1 = @(x) (x-min(x(:)))./(max(x(:))-min(x(:)));

cleanECGScaled = scale0To1(cleanECG);
addNoiseScaled = scale0To1(addNoise);
atrialFibrillationScaled = scale0To1(atrialFibrillation);


%% Generate ECG with undergoing AF

% add signals
alphaAF    = 0.8;
alphaNoise = 0.5;
alphaECG   = 2.5;
simECGwithAF = alphaECG*cleanECGScaled + alphaNoise*addNoiseScaled + alphaAF*atrialFibrillationScaled;

% figure
currentFigName = 'SynthECG';
figTitle = sprintf('%s %s', currentTest, currentFigName);

figure,
subplot(211),
plot(alphaECG*cleanECGScaled),
hold on
plot(alphaNoise*addNoiseScaled, 'g'),
plot(alphaAF*atrialFibrillationScaled, 'r');
legend('synthECG', 'addedNoise', 'synthAF');
title('Synthtesized ECG+AF');
subplot(212),
plot(simECGwithAF)
%set(gcf, 'units', 'normalized', 'outerposition',[0 0 1.5 1.5]);
set(gcf, 'units', 'normalized', 'outerposition',[0 0 1 1]);

figPathName = fullfile(mainFolder, [figTitle, '.png']);
hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
%close(gcf);

% figure,
% subplot(211),
% plot(cleanECG),
% hold on
% plot(addNoise, 'g'),
% plot(10*atrialFibrillation, 'r');
% legend('synthECG', 'addedNoise', 'synthAF');
% title('Synthtesized ECG+AF');
% subplot(212),
% plot(simECGwithAF)
% set(gcf, 'units', 'normalized', 'outerposition',[0 0 1.5 1.5]);
% 
% figPathName = fullfile(mainFolder, [figTitle, '.png']);
% hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
% close(gcf);


%% Prepare data for fastICA

% create mixtures
%inputSignals = [cleanECG; atrialFibrillation; addNoise];
inputSignals   = [alphaECG*cleanECGScaled; alphaAF*atrialFibrillationScaled; alphaNoise*addNoiseScaled];
signalNamesICA = {'synthECG', 'synthAF', 'addedNoise'};

Amatrix  = rand(size(inputSignals,1));
mixedSig = Amatrix*inputSignals;

% figure
currentFigName = 'mixed Signals';
figTitle = sprintf('%s %s', currentTest, currentFigName);

figure,
idx = 1;
subplot(3,1,idx), plot(mixedSig(idx, :));
title('Mixed signals');
idx = 1 + idx;
subplot(3,1,idx), plot(mixedSig(idx, :), 'r')
idx = 1 + idx;
subplot(3,1,idx), plot(mixedSig(idx, :), 'g')


%set(gcf, 'units', 'normalized', 'outerposition',[0 0 1.5 1.5]);
set(gcf, 'units', 'normalized', 'outerposition',[0 0 1 1]);
figPathName = fullfile(mainFolder, [figTitle, '.png']);
hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
%close(gcf);
  

% save data
testData = struct;
testData.simECG             = simECG;
testData.cleanECG           = cleanECG;
testData.addNoise           = addNoise;
testData.atrialFibrillation = atrialFibrillation;
testData.Amatrix            = Amatrix;
testData.mixedSig           = mixedSig;

matFile = [currentTest, '.mat'];
save(fullfile(mainFolder, matFile), 'testData');


%%
% dataIn.signalData = mixedSig;
% dataIn.fs  = 1e3;
% dataOut    = preprocessForICA (dataIn);
% dataForICA = dataOut.signalData;

dataForICA = mixedSig;




%% run fastICA and sort output by excess Kurtosis
[numInputSignalsICA, ~] = size(dataForICA);
[signalICA, A, W] = fastica(dataForICA);

kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);
icaSorted = signalICA(idxKurt, :);

numOutputSignalsICA = numel(kurtosisICA);

if numInputSignalsICA ~= numOutputSignalsICA
  warning('fastICA found %d independent components (%d signals slashed)\n', ...
  numOutputSignalsICA, numInputSignalsICA- numOutputSignalsICA); %#ok<WNTAG>
end




%% prepare data for representation
addOne = @(x) x+1;
addsubplotindex = ...
  @(f,n) setappdata(f, 'subplotindices', [getappdata(f,'subplotindices'), n]);

numPlotColumns = 7;
numPlots       = 5;
offsetCol      = [1,0,1,0,0];
handlesArray   = zeros(1, numPlots*numOutputSignalsICA);

currentIdx     = 1;
currentHandle  = 1;
currentFig     = figure();

for idx = 1:numOutputSignalsICA  
  for idxCol = 1:numPlots      
  handlesArray(currentHandle) = subplot( numOutputSignalsICA, numPlotColumns, currentIdx + [0,offsetCol(idxCol)] );
  currentIdx    = 1 + currentIdx + offsetCol(idxCol);
  currentHandle = addOne(currentHandle);
  end
end

%% survey ICA output and plot

dataIn.fs = 1e3;

normalise = @(x) x./sum(x);

icaFeatures   = struct;
currentHandle = 1;

for idx=1:numOutputSignalsICA

  outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
  
  h0 = handlesArray(currentHandle); 
  plot (h0, dataForICA(idx,:), 'LineWidth', 2);
  title(h0, signalNamesICA{idx});
  axis (h0,'tight');

  currentHandle = addOne(currentHandle);
  h1 = handlesArray(currentHandle);
  
  % get the PSD for the current input signal
  [inputPSD, inputPSDOmega] = getPowerSpectralDensity(dataForICA(idx,:), dataIn.fs, ...
  'windowType', 'rectangular', 'windowSize', length(dataForICA(idx,:)), ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 30);

  inputPSD = normalise(inputPSD);  

  plot(h1, inputPSDOmega, inputPSD, 'LineWidth', 2);
  [~, inputFreqPeakIdx] = max(inputPSD);
  axis(h1, 'tight');
  hold(h1);
  stem(h1, inputPSDOmega(inputFreqPeakIdx), inputPSD(inputFreqPeakIdx), ...
  'LineWidth', 2);
  psdText  = sprintf('PSD F_d = %3.2f Hz', inputPSDOmega(inputFreqPeakIdx));
  title(h1, psdText);   
  set(h1, 'xtick',[]);
  set(h1, 'xticklabel',[]);
  set(h1, 'ytick',[]);
  set(h1, 'yticklabel',[]);   
  
 
  currentHandle = addOne(currentHandle);
  h2 = handlesArray(currentHandle);
  
  currentHandle = addOne(currentHandle);
  h3 = handlesArray(currentHandle);
  currentHandle = addOne(currentHandle);
  h4 = handlesArray(currentHandle); 
  
  [icaFeatures(idx).inputDataStats, icaFeatures(idx).powSpectDens, ...
    icaFeatures(idx).psdOmega, icaFeatures(idx).psdPeaks] = ...
    surveyIndependentComponents( icaSorted(idx, :), outSignalName, dataIn.fs, h2, h3, h4);
  
  set(h2,'ButtonDownFcn',@(~,~) addSubplotIndex(currentFig, idx, h2));
  currentHandle = addOne(currentHandle);
    
end

h = uicontrol('String', 'ICs selected', 'Position', [20 20 100 30], ...
'Callback', 'set(gcbf, ''Name'', sprintf(''IC %d '', getappdata(gcbf,''subplotindices'')))');
% this one doesn't get the figure 'completely' maximised
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
% set(get(handle(currentFig),'JavaFrame'), 'Maximized', 1);
waitfor(currentFig, 'Name');

%end % remove that one

%% Re-construct signals with selected components
%
%  In this step we set to zero all of the sources regarded as VA and
%  reconstruct observations as $hat{x} = A \times s$

selectedICs = getappdata(currentFig, 'subplotindices');

% Get the observations from the sources. 
%   Set to zero all the ICs, and add back the ones picked.
%   Keep the original order as they were displatyed by eKurt order.
selectedSources = zeros(size(signalICA));
selectedSources(idxKurt(selectedICs), :) = icaSorted(selectedICs, :);

%   $hat{x} = A \times s$
dataForICAHat =  A * selectedSources;


% Overlay Xhat obsevations
currentHandle = 1;

for idx=1:numOutputSignalsICA
    
  h0 = handlesArray(currentHandle); 
  hold(h0); %hold current data
  plot(h0, dataForICAHat(idx,:), 'LineWidth', 1, 'Color', 'r');
  
  currentHandle = addOne(currentHandle);
  h1 = handlesArray(currentHandle); 
  [inputHatPSD, inputHatPSDOmega] = getPowerSpectralDensity(dataForICAHat(idx,:), dataIn.fs, ...
  'windowType', 'rectangular', 'windowSize', length(dataForICAHat(idx,:)), ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 30);
  
  % normalise the PDS so we can overlay the Xhat PSD
  inputHatPSD = normalise(inputHatPSD);  

  hold(h1, 'on'); %hold current data
  plot(h1, inputHatPSDOmega, inputHatPSD, 'LineWidth', 1, 'Color', 'r');
  [~, inputHatFreqPeakIdx] = max(inputHatPSD);
  hold(h1, 'on');
  stem(h1, inputHatPSDOmega(inputHatFreqPeakIdx), inputHatPSD(inputHatFreqPeakIdx), ...
  'LineWidth', 1, 'Color', 'r');
  axis(h1, 'tight');

  currentPSDText = get(get(h1, 'title'),'String');% get text back before overwritting
  psdText = sprintf('%s\nIC F_d = %3.2f Hz', currentPSDText, inputHatPSDOmega(inputHatFreqPeakIdx));
  title(h1, psdText);
  
  currentHandle = currentHandle + numPlots-1;  
  
end

% update results
for idx = selectedICs
    icaFeatures(idx).icSelected = 1;
end


% save figure
%set(currentFig, 'units', 'normalized', 'outerposition',[0 0 1.5 1.5]);
set(currentFig, 'units', 'normalized', 'outerposition',[0 0 1 1]);
set(currentFig, 'Name' , figTitle);

% figPathName = fullfile(resultsFolder, [figTitle, '.png']);
% hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
% close(currentFig);
% fprintf('\n');

% 
testResults{currentPoint} = icaFeatures;


%%
idx = 2;
%currentIC = scale0To1(-1*icaSorted(idx, :));
currentIC = scale0To1(icaSorted(idx, :));
figure,
subplot(211)
plot(currentIC)
hold on
plot(atrialFibrillationScaled, 'g')
subplot(212)
plot(abs(currentIC-atrialFibrillationScaled))