function plotCartoAnalysisOverlayCallback(currentHandle)

% Retrieve the analysis data from the current handle and call
% PLOTCARTOANALYSISOVERLAY

analysisData = getappdata(currentHandle, 'cartoAnalysis');  
plotCartoAnalysisOverlay(analysisData);