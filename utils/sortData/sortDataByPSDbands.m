function [inputDataSorted, psdValues, idxPSD] = ...
    sortDataByPSDbands(inputData, samplingFrequency, lowerFreq, upperFreq)

% SORTDATABYPSDBANDS sort data acording to the PDS distribution in a given
% range of frequencies.
%
%

normaliseToUnit   = @(x) x./sum(x);
belowThresholdIdx = @(x, uB) x <= uB;
aboveThresholdIdx = @(x, lB) x >= lB;
withinBoundsIdx   = @(x, lB, uB) aboveThresholdIdx(x, lB) & belowThresholdIdx(x, uB);

[numSignals, numSamples] = size(inputData);
psdWithinBounds = zeros(1, numSignals);


for idx = 1:numSignals
    
    % get the PSD
    [inputPSD, inputPSDOmega] = getPowerSpectralDensity(inputData(idx,:), samplingFrequency, ...
    'windowType', 'rectangular', 'windowSize', numSamples, ...
    'lengthFFT', 8192, 'numberOverlap', [], 'psdCutOff', 30);

    % get the indices for the current band of interest
    currentIdx = withinBoundsIdx(inputPSDOmega, lowerFreq, upperFreq);

    % normalise the PSD to unit and get the normalised energy within bound
    normInputPDS         = normaliseToUnit(inputPSD);
    psdWithinBounds(idx) = sum(normInputPDS(currentIdx));

end


[psdValues, idxPSD] = sort(psdWithinBounds, 'descend');
inputDataSorted     = inputData(idxPSD, :);

% TO-DO: Add here code for svn/gitHub  
% $Revision: $
% $Author: carlosAguilar $
% $Date: $
% Copyright 2015 OGTel LTD.