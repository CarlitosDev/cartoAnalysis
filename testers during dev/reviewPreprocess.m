 idx = 1
 samplingFreq = 1e3;
plotSignalAndSpectrum(dataForICA(idx, :), samplingFreq)
idx = 2
plotSignalAndSpectrum(dataForICA(idx, :), samplingFreq)


idx = 1;
inputECG = dataForICA(idx, :);
idx = 2;
inputEGM = dataForICA(idx, :);
plotEGMandECGTimeAndSpectrum(inputEGM, inputECG, samplingFreq)

idx = 1;
samplingStep = 20;
samplingFreq = 1e3;
inputEGM = dataForICA(idx, :);
[detrendedData, baselineSignal] = removeBaselineWanderingV2(inputEGM, samplingStep);

%plotEGMandECGTimeAndSpectrum(detrendedData, inputEGM, samplingFreq)
plotEGMandECGTimeAndSpectrum(detrendedData, baselineSignal, samplingFreq)


figureLeft, plot(padInputData)
hold on
plot(filteredInputData, 'r')

figureRight, 
plot(filteredInputData), 
hold on, plot(padBaselineSignal, 'r')
hold on, stem(xControlPoints, yControlPoints)
