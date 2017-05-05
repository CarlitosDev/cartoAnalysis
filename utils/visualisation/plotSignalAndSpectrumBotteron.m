function plotSignalAndSpectrumBotteron(inputSignal, samplingFreq)

% do it right at some time later

% create a figure and plot signal in time domain

 figure(),
  
  subplot(211),
  plot(inputSignal, 'Color', [1 0 0], 'LineWidth', 1.5);
  
  % get the PSD for the current input signal
  [inputPSD, inputPSDOmega] = getPowerSpectralDensityBotteron(inputSignal, samplingFreq, ...
  'windowType', 'rectangular', 'windowSize', length(inputSignal), ...
  'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 60);
  
  subplot(212),
  plot(inputPSDOmega, inputPSD, 'LineWidth', 2);


