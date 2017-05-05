samplingFreq = 1e3;
inputSignal = cartoPointAnalysis{1, 1}.dataForICA;
%botteronPreprocessed = preprocessBotteronSmith(inputSignal, samplingFreq);

lengthSignal = numel(inputSignal);

[inputHatPSD, inputHatPSDOmega] = ...
    getPowerSpectralDensityBotteron(inputSignal, samplingFreq, ...
    'windowType'   , 'rectangular', ...
    'windowSize'   , lengthSignal/4, ...
    'lengthFFT'    , 8192, ...
    'numberOverlap', [], ...
    'psdCutOff'    , 60);

figure,
plot(inputHatPSDOmega, inputHatPSD);

% figure,
% plot(inputSignal, 'LineWidth', 1.5, 'Color', [1 0 0]);
% hold on
% plot(dataSignalProcWhite, 'LineWidth', 1.5, 'Color', [0 1 0]);
