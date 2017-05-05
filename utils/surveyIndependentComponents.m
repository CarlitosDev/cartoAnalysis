function [inputDataStats ,powSpectDens, psdOmega, psdPeaks] = ...
    surveyIndependentComponents( inputData, inputDataName, inputDataFs, h1, h2, h3)
  

%%   SURVEYINDEPENDENTCOMPONENTS get information about an input IC.

plotSecondHarmonic  = false;
useWelchPeriodogram = false;

normalise = @(x) x./sum(x);


fprintf(' Stats for signal %s', inputDataName) ;
inputDataStats = Utils_basic_stats( inputData );

% get PSD
getQuarterPow2    = @(x) 2^(nextpow2(length(x))-2);
hammingWindowSize = getQuarterPow2(inputData);

if useWelchPeriodogram
    % Hamming window
    [powSpectDens, psdOmega] = getPowerSpectralDensity(inputData, inputDataFs, ...
      'windowSize', hammingWindowSize, 'lengthFFT', 8192,  ...
      'numberOverlap', [], 'psdCutOff', 0.9); %#ok<UNRCH>
else
    % Rectangular window
    [powSpectDens, psdOmega] = getPowerSpectralDensity(inputData, inputDataFs, ...
      'windowType', 'rectangular', ...
      'windowSize', hammingWindowSize, 'lengthFFT', 8192,  ...
      'numberOverlap', [], 'psdCutOff', 20);
end

% plot
if ishandle(h1)
plot (h1, inputData, 'LineWidth', 2, 'Color', [0 0 1]);
title(h1, inputDataName);
axis (h1,'tight')
set  (h1,'xtick',[])
set  (h1,'xticklabel',[])
end

% PSD. Get 2 peaks
if ishandle(h2)
    plot(h2, psdOmega, powSpectDens, 'LineWidth', 2, 'Color', [0 0 1]);
    axis(h2,'tight')
    %set(h2 , 'xtick', [])
    %set(h2 , 'xticklabel', [])   
end    

    [~, idxSorted] = sort(powSpectDens, 'descend');

    indexLTthreshold  = @(x,y,threshold) find(abs(x-y)>threshold, 1, 'first');
    deltaFreq         = 1.5; % Hz
    psdOmegaSorted    = psdOmega(idxSorted);

    secondHarmonicIdx = indexLTthreshold(psdOmegaSorted(2:end), psdOmegaSorted(1), deltaFreq) + 1;
    secondHarmonic    = psdOmega(idxSorted(secondHarmonicIdx));

    secondHarmonicIdx = [];

    % thirdHarmonicIdx  = indexLTthreshold(psdOmegaSorted(3:end), secondHarmonic, deltaFreq) + 2;
    % harmonicsIdx      = [idxSorted(1), idxSorted(secondHarmonicIdx), idxSorted(thirdHarmonicIdx)];

    harmonicsIdx      = [idxSorted(1), idxSorted(secondHarmonicIdx)];
    if ishandle(h2)        
        hold(h2, 'on'); 
        colours = ['r', 'g'];
        if ~isempty(secondHarmonicIdx) && plotSecondHarmonic
            for i=1:2
                stem(h2, psdOmega(harmonicsIdx(i)), powSpectDens(harmonicsIdx(i)), colours(i));
            end
        else
            i = 1;
            stem(h2, psdOmega(harmonicsIdx(i)), powSpectDens(harmonicsIdx(i)), colours(i));
        end
    end

    psdPeaks = psdOmega(harmonicsIdx);
    psdText  = sprintf('F_d = %3.2f Hz', psdPeaks(1));
    
    if ishandle(h2)      
        title(h2,psdText);
        axis(h2,'tight');
    end
% histogram

statsText = sprintf('\\mu = %3.2f, \\sigma = %3.2f', ...
            inputDataStats.Mu, inputDataStats.Std);

numPoints = 25;
xAxis    = linspace(inputDataStats.minVal, inputDataStats.maxVal, numPoints);
normDist = normpdf(xAxis ,inputDataStats.Mu, inputDataStats.Std);
normDist = normDist/sum(normDist);

[freq, freqBin] = hist(inputData, numPoints);
freqNorm = freq/sum(freq);

if ishandle(h3)
    bar(h3, freqBin, freqNorm)
    hold(h3, 'on');
    plot(h3, xAxis, normDist,'r')
    axis(h3,'tight')
    title(h3, statsText);
end