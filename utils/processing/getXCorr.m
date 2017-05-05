function currentData = getXCorr(bulkData, maximumLag)

numSignals  = size(bulkData, 1);
optcor      = 'coeff';
currentData = zeros(numSignals, 1+2*maximumLag);

for idx=1:numSignals
  aux = bulkData(idx, :);
  currentData(idx, :) = xcorr(aux,aux,maximumLag,optcor);
end