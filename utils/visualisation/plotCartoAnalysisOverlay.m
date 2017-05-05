function plotCartoAnalysisOverlay(analysisData)

% PLOTCARTOANALYSISOVERLAY
% Plot the results of the carto analysis. Overlay the EGM with the
% reference EGC on the top figure and show both frequency peak and
% oraanisation index in the bottom one.
  

%% scale input signals
    scaleFrom0To1 = @(x) (x-min(x))/(max(x)-min(x));
    % time signals
    scaledECG     = scaleFrom0To1(analysisData.referenceSignal);
    scaledEGM     = scaleFrom0To1(analysisData.dataForICA);
    centeredECG   = scaledECG - mean(scaledECG, 2);
    centeredEGM   = scaledEGM - mean(scaledEGM, 2);
    % also scale frequency
    inputPSD      = scaleFrom0To1(analysisData.inputPSD);
    inputHatPSD   = scaleFrom0To1(analysisData.inputHatPSD);
    
    frequencyPeak = analysisData.inputHatPSDOmega(analysisData.inputHatFreqPeakIdx);
    

%% plot    

    % colours
    silverColour    = [192,192,192]./255;
    lightCyanColour = [154,255,255]./255;
    orangeColour    = [255 165 0]  ./255;
    

    currentFigure = figure();
    
    subplot(211),
    figATitle = sprintf('Analysis of point %d', analysisData.pointId);
    plot(centeredECG, 'Color', silverColour, 'LineWidth', 0.8);
    hold on,
    plot(centeredEGM, 'Color', [1 0 0], 'LineWidth', 1.5);
    
    legend({'Reference V_1', 'EGM'});
    title(gca, figATitle);

    
    subplot(212),    
    
    plot(analysisData.inputPSDOmega, inputPSD, ...
        'Color', silverColour, 'LineWidth', 1.0, ...
        'LineStyle', '-');
    hold on,
    plot(analysisData.inputHatPSDOmega, inputHatPSD, ...
        'LineWidth', 2);
    
    hold on,
    
    % main frequency
    stem(analysisData.inputHatPSDOmega(analysisData.inputHatFreqPeakIdx), ...
      inputHatPSD(analysisData.inputHatFreqPeakIdx), ... 
      'LineWidth', 2, 'Color', [0 0 0]);
  
    searchWindow = 0.5*analysisData.searchWindow;
    
    leftSideIdx  = analysisData.inputHatPSDOmega >= frequencyPeak-searchWindow;
    rightSideIdx = analysisData.inputHatPSDOmega <= frequencyPeak+searchWindow;
    searchIdx    = leftSideIdx & rightSideIdx;
    
    hold on,
    
    area(analysisData.inputHatPSDOmega(searchIdx), inputHatPSD(searchIdx), ...
        'FaceColor', orangeColour);
    
    legend({'Original', 'Botteron', 'Peak'});
    

    figBTitle = sprintf('F_d %3.2f Hz [%3.2f, %3.2f] - OI  %3.2f', analysisData.inputHatFreqPeak, ...
        frequencyPeak-searchWindow, frequencyPeak+searchWindow, ...
        analysisData.organisationIndex);
    title(figBTitle);

    % Align to the right
    alignFigure(currentFigure, 'right')