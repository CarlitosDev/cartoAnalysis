function plotCartoAnalysisCallback(currentHandle)


% Retrieve the data from the analysis
analysisData = getappdata(currentHandle, 'cartoAnalysis');


currentFigure = figure();

subplot(211),
plot(analysisData.dataForICA, 'Color', [1 0 0], 'LineWidth', 1.5);

subplot(212),
plot(analysisData.inputHatPSDOmega, analysisData.inputHatPSD, 'LineWidth', 2);
hold on,
stem(analysisData.inputHatPSDOmega(analysisData.inputHatFreqPeakIdx), ...
  analysisData.inputHatPSD(analysisData.inputHatFreqPeakIdx)        , ... 
  'LineWidth', 2, 'Color', 'red');

figTitle = sprintf('Point %d peaking at %3.2f Hz', analysisData.pointId, analysisData.inputHatFreqPeak);
title(figTitle);

% Align to the right
alignFigure(currentFigure, 'right')