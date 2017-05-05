function plotDescriptorOverPoints(cathererTable, pointDescriptor, figTitle)

%% PLOTDESCRIPTOROVERPOINTS plot a CartoXP data analysis

%%

numPoints = height(cathererTable);

% arrange catherer data
xData          = cathererTable{:, 'x'};
yData          = cathererTable{:, 'y'};
zData          = cathererTable{:, 'z'};

%% Heatmap colours
% Map the values of the descriptor into the 64 values of a colormap
maxDesc        = max(pointDescriptor);
minDesc        = min(pointDescriptor);
numBins        = 64;
binRanges      = linspace(minDesc, maxDesc, numBins);
[~, colourIdx] = histc(pointDescriptor, binRanges);

%% Plot
cFigure = figureRight();
cMap    = colormap(cFigure, hot);

descriptorColours = cMap(colourIdx, :);


% Matlab doesn't let me pass a Nx3 vector, so let's iterate though
% the points
for iPoint = 1:numPoints
  plot3(xData(iPoint), yData(iPoint), zData(iPoint), ...
    'o', 'MarkerSize'  , 8, ...
    'MarkerEdgeColor', 'none', ... 
    'MarkerFaceColor'  , descriptorColours(iPoint, :));
  hold on;
  if mod(iPoint, 500) == 0
    fprintf('Plotting %d/%d...\n', iPoint, numPoints);
  end
end

cAx = gca();
caxis(cAx, [minDesc, maxDesc])
colorbar(cAx);

title(figTitle);