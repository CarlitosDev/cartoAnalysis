function [inputDataStats,powSpectDens, psdOmega] = ...
    surveyDataForICA( inputData, inputDataName, inputDataFs)

fprintf(' Stats for signal %s', inputDataName) ;
inputDataStats = Utils_basic_stats( inputData );

% get PSD
getQuarterPow2    = @(x) 2^(nextpow2(length(x))-2);
hammingWindowSize = getQuarterPow2(inputData);

[powSpectDens, psdOmega] = getPowerSpectralDensity(inputData, inputDataFs, ...
  'hammingWindow', hammingWindowSize, 'lengthFFT', 8192,  ...
  'numberOverlap', [], 'psdCutOff', 0.98);


% plot
figure,
subplot(131), plot(inputData);
ylabel(inputDataName);
statsText = ...
    sprintf('\\mu = %3.2f, \\sigma = %3.2f, median = %3.2f. Range (%3.2f, %3.2f)', ...
    inputDataStats.Mu    , inputDataStats.Std, ...
    inputDataStats.median, ...
    inputDataStats.minVal, inputDataStats.maxVal);
xlabel(statsText);


% PSD. Get 2 peaks
subplot(132), 
plot(psdOmega, powSpectDens)
[~, idxSorted] = sort(powSpectDens, 'descend');

indexLTthreshold  = @(x,y,threshold) find(abs(x-y)>threshold, 1, 'first');
deltaFreq         = 1; % Hz
psdOmegaSorted    = psdOmega(idxSorted);

secondHarmonicIdx = indexLTthreshold(psdOmegaSorted(2:end), psdOmegaSorted(1), deltaFreq) + 1;
secondHarmonic    = psdOmega(idxSorted(secondHarmonicIdx));

% thirdHarmonicIdx  = indexLTthreshold(psdOmegaSorted(3:end), secondHarmonic, deltaFreq) + 2;
% harmonicsIdx      = [idxSorted(1), idxSorted(secondHarmonicIdx), idxSorted(thirdHarmonicIdx)];

harmonicsIdx      = [idxSorted(1), idxSorted(secondHarmonicIdx)];
hold on,
colours = ['r', 'g'];
for i=1:2
    stem(psdOmega(harmonicsIdx(i)), powSpectDens(harmonicsIdx(i)), colours(i));
end
psdText = sprintf('max F_d = %3.2f Hz, F_{d2} = %3.2f', ...
    psdOmega(harmonicsIdx(1)), psdOmega(harmonicsIdx(2)));
xlabel(psdText);

% histogram
numPoints = 25;
xAxis    = linspace(inputDataStats.minVal, inputDataStats.maxVal, numPoints);
normDist = normpdf(xAxis ,inputDataStats.Mu, inputDataStats.Std);
normDist = normDist/sum(normDist);

subplot(133), 
[freq, freqBin] = hist(inputData, numPoints);
freqNorm = freq/sum(freq);
bar(freqBin, freqNorm)
hold on,
plot(xAxis, normDist,'r')