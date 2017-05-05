
numPoints = numel(cartoPointAnalysis);

  inputHatFreqPeak  = zeros(1, numPoints);
  organisationIndex = zeros(1, numPoints);
  maxVolt           = zeros(1, numPoints);
  minVolt           = zeros(1, numPoints);
  meanVolt          = zeros(1, numPoints);

for idxCarto = 1:numPoints
  
  cPoint = cartoPointAnalysis{idxCarto};

  inputHatFreqPeak(idxCarto)  = cPoint.inputHatFreqPeak;
  organisationIndex(idxCarto) = cPoint.organisationIndex;
  maxVolt(idxCarto)  = max(cPoint.dataForICA);
  minVolt(idxCarto)  = min(cPoint.dataForICA);
  meanVolt(idxCarto) = mean(cPoint.dataForICA);
end

%%
cData = [inputHatFreqPeak; organisationIndex; ...
   meanVolt];

pointNames = regexp(sprintf('P_%d\n', 1:numPoints), '\n', 'split');
pointNames(end) = [];

yValues = {'inputHatFreqPeak', 'organisationIndex',  ...
  'meanVolt'};
figure,
plot(1:numPoints, cData, 'LineWidth', 2);
legend(yValues);

%%




h = heatmap(pointNames, yValues, cData);


%%
numPoints = 200
pointNames = regexp(sprintf('P_%d\n', 1:numPoints), '\n', 'split');
pointNames(end) = [];

yValues = {'inputHatFreqPeak', 'organisationIndex', 'maxVolt',  ...
  'minVolt','meanVolt'};

cData = [inputHatFreqPeak; organisationIndex; maxVolt; ...
  minVolt; meanVolt];

cData = cData(:, 1:numPoints);

h = heatmap(pointNames, yValues, cData);

h.Title = 'Parameters';
h.XLabel = 'Points';
h.YLabel = 'Descriptors';

%%

figure,
plot(1:numPoints, cData);
legend(yValues);

%%


