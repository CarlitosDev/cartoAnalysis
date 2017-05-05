cathTable    = struct2table(cathererData);
idxElectrode = cathTable.pentaRayPoint == targetElectrode;
cathTable    = cathTable(idxElectrode, :);
cathTable.pointsId(:) = 1:height(cathTable);



idxSignal = find(contains(signalNames, 'x20A_1_31_'));

figure,
plot(tempData(idxPoint).data{idxSignal})

hold on,
idxSignal = find(contains(signalNames, 'x20A_2_32_'));

plot(tempData(idxPoint).data{idxSignal}, 'r')

hold on,
idxSignal = find(contains(signalNames, 'x20A_3_33_'));
plot(tempData(idxPoint).data{idxSignal}, 'g')


hold on,
idxSignal = find(contains(signalNames,  'x20A_5_35_'));
plot(tempData(idxPoint).data{idxSignal}, '-b')

%%
currentString = 'x20A_1_2_31_'
idxA = strfind(currentString, '_')
currentString(idxA(end-1):end) = []

%%

plot3(ePos.X, ePos.Y, ePos.Z, 'o')
for i=1:height(ePos)

 text(ePos.X(i)+0.15, ePos.Y(i), ePos.Z(i), num2str(ePos.electrodeNumbers(i)))
 hold on
end

%%
cathTable.Properties.VariableNames'

listOfPoints = unique(cathTable{:, 'pentaRayPoint'});
idxSubSet    = cathTable{:, 'pentaRayPoint'} == listOfPoints(50);

subSet = cathTable(idxSubSet, :);
idxMon = strcmpi(subSet{:, 'type'}, 'BIP');

subSet = subSet(idxMon, :);

hFig = figureMax();
lAx  = subplot(2,2,[1 3]);
hold(lAx, 'on');
plot3(lAx, subSet{:, 'x'}, ...
  subSet{:, 'y'}, subSet{:, 'z'}, 'o');
view(-10, 30);
title(lAx, 'Catheter position');

rAx1 = subplot(2,2,2);
hold(rAx1, 'on');
rAx2 = subplot(2,2,4);
hold(rAx2, 'off');

numSignals = height(subSet);

currentShift = 0;

for idx=1:numSignals
  
  text(lAx, 0.2+subSet{idx, 'x'}, ...
    subSet{idx, 'y'}, subSet{idx, 'z'}, ...
    subSet{idx, 'pentaRayElectrode'});
  
  plot(rAx1 , currentShift+subSet{idx, 'signal'});
  currentShift = currentShift + max(subSet{idx, 'signal'});
  plot(rAx2 , subSet{idx, 'signal'});
  title(rAx2, subSet{idx, 'pentaRayElectrode'});
  pause(2)
  
end



%%
cathTable.Properties.VariableNames'

listOfPoints = unique(cathTable{:, 'pentaRayPoint'});
idxSubSet    = cathTable{:, 'pentaRayPoint'} == listOfPoints(20);

subSet = cathTable(idxSubSet, :);
idxMon = strcmpi(subSet{:, 'type'}, 'BIP');

subSet = subSet(idxMon, :);

hFig = figureMax();
lAx  = subplot(1,2,1);
hold(lAx, 'on');
plot3(lAx, subSet{:, 'x'}, ...
  subSet{:, 'y'}, subSet{:, 'z'}, 'o');
view(-10, 30);
title(lAx, 'Catheter position');

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
    10-currentShift+mean(currentSignal(end-200:end-100)), ...
    subSet{idx, 'pentaRayElectrode'});
  
  currentShift = currentShift + max(subSet{idx, 'signal'});
  
end
