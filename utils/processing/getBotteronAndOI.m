
function [botFrequencies, organisationIdx] = ...
  getBotteronAndOI(inputData, samplingFreq, oiBandWidth)

nR = size(inputData, 1);

botFrequencies  = zeros(1, nR);
organisationIdx = zeros(1, nR);

for idx=1:nR
  
  currentSignal = inputData(idx, :);
  [botFreq,~,~,~,periodogramInput,fz] = df_Ng(currentSignal, samplingFreq);
  oi = organizationIndex(periodogramInput,fz,botFreq,oiBandWidth);
  
  botFrequencies(idx)  = botFreq;
  organisationIdx(idx) = oi;
  if mod(idx, 1000) == 0
    fprintf('%d/%d processed...\n', idx, nR)
  end
end