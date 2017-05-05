function [currentData, botFreqs] = getBotteronSignals(bulkData, samplingFreq)

[numSignals, numSamples]  = size(bulkData);
currentData = zeros(numSignals, numSamples);
botFreqs    = zeros(numSignals, 1);

for idx=1:numSignals
  aux = bulkData(idx, :);
  [botFreq,~,~,bottSignal] = df_Ng(aux, samplingFreq);
  currentData(idx, :) = bottSignal;
  botFreqs(idx)       = botFreq;
end