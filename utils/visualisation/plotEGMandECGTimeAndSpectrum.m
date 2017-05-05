function plotEGMandECGTimeAndSpectrum(inputEGM, inputECG, samplingFreq)


%% scale input signals

scaleFrom0To1 = @(x) (x-min(x))/(max(x)-min(x));
scaledECG     = scaleFrom0To1(inputECG);
scaledEGM     = scaleFrom0To1(inputEGM);

centerdECG    = scaledECG - mean(scaledECG, 2);
centerdEGM    = scaledEGM - mean(scaledEGM, 2);

%%
% light cyan
ecgColour = [154,255,255]./255;

figure(),
subplot(211),
plot(centerdECG, 'Color', ecgColour, 'LineWidth', 0.8);
hold on,
plot(centerdEGM, 'Color', [1 0 0], 'LineWidth', 1.8);



% get the PSD for the signals
[egmPSD, egmPSDOmega] = getPowerSpectralDensity(centerdEGM, samplingFreq, ...
'windowType', 'rectangular', 'windowSize', length(centerdEGM), ...
'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 60);

[ecgPSD, ecgPSDOmega] = getPowerSpectralDensity(centerdECG, samplingFreq, ...
'windowType', 'rectangular', 'windowSize', length(centerdECG), ...
'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 60);



subplot(212),
plot(egmPSDOmega, egmPSD, 'Color', [1 0 0], 'LineWidth', 1.5);
hold on,
plot(ecgPSDOmega, ecgPSD, 'Color', ecgColour, 'LineWidth', 1.0);
