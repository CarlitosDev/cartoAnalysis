function plotSignalAndSpectrum(inputSignal, samplingFreq)

% do it right at some time later


% dfSignal = 100; %Hz
% fs = 1000;
% t = 0:1/fs:5-1/fs;
% inputSignal = cos(2*pi*dfSignal*t)+randn(size(t));
% 
% plotSignalAndSpectrum(inputSignal, fs)

% create a figure and plot signal in time domain

 figure(),
  
  subplot(211),
  plot(inputSignal, 'Color', [1 0 0], 'LineWidth', 1.5);
  drawnow;
  
  % get the PSD for the current input signal
  [inputPSD, inputPSDOmega] = getPowerSpectralDensity(inputSignal, samplingFreq, ...
  'windowType', 'rectangular', 'windowSize', length(inputSignal), ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 180);
  
  subplot(212),
  plot(inputPSDOmega, inputPSD, 'LineWidth', 2);


