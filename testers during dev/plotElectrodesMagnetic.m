load('ECG.mat')

ePos = ECG(1).electrodesPos;

figure,
plot3(ePos.X, ePos.Y, ePos.Z, 'o')
text(ePos.X, ePos.Y, 0.25+ePos.Z, num2str([1:numel(ePos.X)]'))

hold on
ePos = ECG(100).electrodesPos;
plot3(ePos.X, ePos.Y, ePos.Z, 'ro')
text(ePos.X, ePos.Y, 0.25+ePos.Z, num2str(ePos.electrodeNumbers))


%
   electroNames = {'NA', ...
    'NA', ...
    'x20A_1_31_', ...
    'x20A_2_32_', ...
    'x20A_3_33_', ...
    'x20A_4_34_', ...
    'x20A_5_35_', ...
    'x20A_6_36_', ...
    'x20A_7_37_', ...
    'x20A_8_38_', ...
    'x20A_9_39_', ...
    'x20A_10_40_', ...
    'x20A_11_41_', ...
    'x20A_12_42_', ...
    'x20A_13_43_', ...
    'x20A_14_44_', ...
    'x20A_15_45_', ...
    'x20A_16_46_', ...
    'x20A_17_47_', ...
    'x20A_18_48_', ...
    'x20A_19_49_', ...
    'x20A_20_50_'};

  for i=1:22
  electroNames{i} = sprintf('Point %d (%s)', ePos.electrodeNumbers(i), electroNames{i});
  end
  
  figureRight
  ePos = ECG(100).electrodesPos;
plot3(ePos.X, ePos.Y, ePos.Z, 'ro')
  text(ePos.X, ePos.Y, 0.25+ePos.Z, electroNames, 'Interpreter', 'none')

  
%% Plot the data in cartoData against the MAGNETIC catheter

[data, pointsInfo, meshData, cartoData] = loadRamonyCajalData(fullPath2CartoFile);

figureRight,

plot3(cartoData.x, cartoData.y, cartoData.z, ...
  'o', 'MarkerSize' , 4, ...
  'MarkerFaceColor' , [0 1 0], ...
  'MarkerEdgeColor' , 'none');

numAcqPoints = numel(ECG);
drawnow

for idxPoint=1:numAcqPoints
  hold on,
  ePos = ECG(idxPoint).electrodesPos;
  plot3(ePos.X, ePos.Y, ePos.Z, 'ro');
  pause(0.05);
end

%%

idxElectrode = 1;

  ePos = ECG(idxElectrode).electrodesPos;

distPoint2PointSet = ...
    @(pointSet, point) sqrt(sum(bsxfun(@minus, pointSet, point).^2, 2));
  
currentPointSet = [cartoData.x; cartoData.y; cartoData.z]';  
currentPoint    = [ePos.X, ePos.Y, ePos.Z];


[numEletrodes, ~] = size(currentPoint);
for idxPoint = 1:numEletrodes
	[pointFacetDist, pointIdx] = ...
      min(distPoint2PointSet(currentPointSet, currentPoint(idxPoint, :)));
   if  pointFacetDist < 1e-2
     fprintf('Catherer point %d is carto point %d\n', idxPoint, pointIdx);
   end
    
end

    
%     plot3(currentPoint(idxPoint, 1), currentPoint(idxPoint, 2), currentPoint(idxPoint, 3), 'ro');