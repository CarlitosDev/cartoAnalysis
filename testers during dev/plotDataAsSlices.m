function plotDataAsSlices(data2plot, xAxisPosition, xAxisTickName)



% 
% data2plot     = rand(4, 1000);
% xAxisPosition = randi(30, [1,4]);
% xAxisTickName = {'P1', 'P2', 'P3', 'P4'};


[numSignals, numSamples] = size(data2plot);

yAxis = repmat(1:numSamples, [numSignals,1]);

%%

h = figure();
hAx = gca();
for idx = 1:numSignals
  plot3(hAx, repmat(xAxisPosition(idx), [1, numSamples]), yAxis(idx, :), data2plot(idx, :));
  hold(hAx, 'on')
end

currentXTick       = get(hAx, 'XTick');
currentXTickLabels = get(hAx, 'XTickLabels');

set(hAx, 'XTickLabels', xAxisTickName);