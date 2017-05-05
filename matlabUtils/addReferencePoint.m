function addReferencePoint(currentFig, currentPointIdx)

currentIndexes = getappdata(currentFig, 'pointIndices');
currentIndexes = unique([currentIndexes, currentPointIdx]);

setappdata(currentFig, 'pointIndices', currentIndexes);

% TO-Do: Work out the index within the Line vector
% hLine = findall(currentFig,'type', 'line');
% currentPlot = hLine(currentPointIdx);
% currentPlot.MarkerSize      = 20;
% currentPlot.MarkerFaceColor = [1 1 0];
% drawnow;

currentTitle = sprintf('%s %d %s', 'Point', currentPointIdx, ...
  'Selected. Press continue');
title(currentTitle);

fprintf('Adding pointIdx %d\n', currentPointIdx);