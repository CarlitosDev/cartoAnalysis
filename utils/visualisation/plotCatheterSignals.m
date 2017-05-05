function plotCatheterSignals(cathTable, currentPoint, currentType)

% PLOTCATHETERSIGNALS plot Pentaray spatial and temporal.
%
% cathTable: A table containing the catheter points. Can be generated with 'loadRamonYCajalData.m'
% 
% currentPoint: A valid pointId from the Pentaray list
% 
% currentType: Either 'MON' or 'BIP'


idxSubSet = cathTable{:, 'pentaRayPoint'} == currentPoint;
subSet    = cathTable(idxSubSet, :);
idxMon    = strcmpi(subSet{:, 'type'}, currentType);
subSet    = subSet(idxMon, :);


hFig = figureMax();
lAx  = subplot(1,2,1);
hold(lAx, 'on');
plot3(lAx, subSet{:, 'x'}, ...
  subSet{:, 'y'}, subSet{:, 'z'}, 'o');
view(-10, 30);
title(lAx, sprintf('%d - %s catheter placement',currentPoint, currentType));

rAx1 = subplot(1,2,2);
hold(rAx1, 'on');

numSignals   = height(subSet);
currentShift = 0;

for idx=1:numSignals
  
  text(lAx, 0.2+subSet{idx, 'x'}, ...
    subSet{idx, 'y'}, subSet{idx, 'z'}, ...
    subSet{idx, 'pentaRayElectrode'});
  currentSignal = subSet{idx, 'signal'};
  plot(rAx1, currentShift+currentSignal);

  
	text(rAx1, length(currentSignal)-300, ...
    -10+currentShift+mean(currentSignal(end-200:end-100)), ...
    subSet{idx, 'pentaRayElectrode'});
  
  currentShift = currentShift + max(subSet{idx, 'signal'});
  
end