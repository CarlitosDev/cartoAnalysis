%% run fastICA and sort output by excess Kurtosis
[signalICA, A, W] = fastica(dataForICA);


a = dataForICA(1, :);


b = W*dataForICA;

figure,
plot(signalICA(1, :), 'r')
hold on
plot(b(1, :))

c = W(:, 1)*dataForICA(1, :);

d = A*W


% recover sources
dataForICAHat = A*signalICA;

figure,
plot(dataForICA(1, :), 'r')
hold on
plot(dataForICAHat(1, 1:1000))


e = signalICA(1,:);
dataForICAHat2 = e*A;

size(e)
size(A)

%
signalICA(1,1)

s11 = sum(W(1,:).*dataForICA(:,1)')
s1N = sum(W(1,:).*dataForICA(:,end)')


%%%%%

[signalICA, A, W] = fastica(dataForICA);
% recover sources
dataForICAHat = A*signalICA;

figure,
plot(dataForICA(1, :), 'r')
hold on
plot(dataForICAHat(1, 1:1000))


%%%%% WRONG order

[signalICA, A, W] = fastica(dataForICA);

kurtosisICA = kurtosis(signalICA, 1, 2) - 3;

[eKurtValues, idxKurt] = sort(kurtosisICA);
icaSorted = signalICA(idxKurt, :);

% recover sources
dataForICAHat = zeros(size(signalICA));
dataForICAHat(idxKurt, :) = icaSorted;
dataForICAHat = A*dataForICAHat;


figure,
plot(dataForICA(1, :), 'r')
hold on
plot(dataForICAHat(1, 1:1000))


%% %%% WRONG order

[signalICA, A, W] = fastica(dataForICA);

% prepare data for representation
addOne = @(x) x+1;
addsubplotindex = ...
  @(f,n) setappdata(f, 'subplotindices', [getappdata(f,'subplotindices'), n]);

numPlotColumns = 6;
numPlots       = 4;
offsetCol      = [1,1,0,0];
handlesArray   = zeros(1, numPlots*numOutputSignalsICA);

currentIdx     = 1;
currentHandle  = 1;
currentFig     = figure();

for idx = 1:numOutputSignalsICA  
  for idxCol = 1:4       
  handlesArray(currentHandle) = subplot( numOutputSignalsICA, numPlotColumns, currentIdx + [0,offsetCol(idxCol)] );
  currentIdx    = 1 + currentIdx + offsetCol(idxCol);
  currentHandle = addOne(currentHandle);
  end
end

% survey ICA output and plot


icaFeatures   = struct;
currentHandle = 1;

for idx=1:numOutputSignalsICA

  outSignalName = sprintf('ICA #%d k = %3.2f', idx, eKurtValues(idx));
  
  h0 = handlesArray(currentHandle); 
  plot (h0, dataForICA(idx,:))
  title(h0, signalNamesICA{idx});
  axis (h0,'tight')
  
  currentHandle = addOne(currentHandle);
  h1 = handlesArray(currentHandle);
  
  currentHandle = addOne(currentHandle);
  h2 = handlesArray(currentHandle);
  currentHandle = addOne(currentHandle);
  h3 = handlesArray(currentHandle); 
  
  [icaFeatures(idx).inputDataStats, icaFeatures(idx).powSpectDens, ...
    icaFeatures(idx).psdOmega, icaFeatures(idx).psdPeaks] = ...
    surveyIndependentComponents( icaSorted(idx, :), outSignalName, dataIn.fs, h1, h2, h3);
  
  set(h1,'ButtonDownFcn',@(~,~) addSubplotIndex(currentFig, idx, h1));
  currentHandle = addOne(currentHandle);
    
end


h = uicontrol('String', 'ICs selected', 'Position', [20 20 100 30], ...
'Callback', 'set(gcbf, ''Name'', sprintf(''IC %d '', getappdata(gcbf,''subplotindices'')))');
waitfor(currentFig, 'Name');


selectedICs = getappdata(currentFig, 'subplotindices');


% recover sources
dataForICAHat = zeros(size(signalICA));
dataForICAHat(idxKurt(selectedICs), :) = icaSorted;


dataForICAHat = A*dataForICAHat;


figure,
plot(dataForICA(1, :), 'r')
hold on
plot(dataForICAHat(1, 1:1000))

