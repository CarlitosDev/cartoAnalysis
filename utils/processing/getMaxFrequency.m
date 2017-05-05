
function maxFrequencies = getMaxFrequency(inputData, samplingFreq, fCutOff)

nR = size(inputData, 1);

maxFrequencies = zeros(1, nR);
numSamples     = numel(inputData(1, :));

for idx=1:nR
  
  currentSignal = inputData(idx, :);
  
  [currentFreqValues, currentFreqAxes] = ...
    getPowerSpectralDensity(currentSignal, samplingFreq, ...
    'windowType', 'rectangular', ...
    'windowSize', numSamples, ...
    'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', fCutOff);
  
  [~, idxPeak] = max(currentFreqValues);
  
  maxFrequencies(idx)  = currentFreqAxes(idxPeak);
  if mod(idx, 100) == 0
    fprintf('%d/%d processed...\n', idx, nR)
  end
end