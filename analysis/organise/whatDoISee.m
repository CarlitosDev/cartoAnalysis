%% Set up


testType = 'simBotteronSmith';

% Keep the simulation off GIT
%baseUngittedFolder = 'D:\Work\MATLABPromotionsModel\Matlab PhD';
baseUngittedFolder = ...
  '/Users/carlosAguilar/Documents/Matlab PhD/Data (ungitted)';


%% Load simulation
fullPath2SimulationFile = fullfile(baseUngittedFolder,  ...
  'SimulacionValencia', 'AF_AI.mat');
[data, pointsInfo, meshData, cartoData] = loadSimulationData(fullPath2SimulationFile);


%%




% 
% fullPath2SimulationFile = fullfile(baseUngittedFolder,  ...
%   'SimulacionValencia', 'analysis17Feb.mat');
% 
% 
% load(fullPath2SimulationFile);

%% Plt a time slice - decimated version
% timeFrame = 3333;
timeFrame = 10;

pointDescriptor = zeros(size(cartoData.pointsId));

for idx=1:numel(cartoData.pointsId)
  pointDescriptor(idx) = data(idx).signal{1,1}(timeFrame);
end

plotMeshAndEGMvalues(meshData, cartoData, pointDescriptor);

%% Compare against all points
%

meshVertex = meshData.meshVertex;
egmData    = cartoData.egmData;
facets     = meshData.faces;

x = 35.0*double(meshVertex(:, 1));
y = 35.0*double(meshVertex(:, 2));
z = 35.0*double(meshVertex(:, 3));

facetColours = egmData(timeFrame, :);

figureLeft,
trisurf(facets, x,y,z, ...
facetColours, ...
'facealpha' , 0.35 )

%% Quick lookup
idx = 100;
figureRight,
plot(data(idx).signal{1,1}, 'r')

idx = 20;
figureLeft,
plot(data(idx).signal{1,1}, 'b')

%%
samplingFreq = 500;

idx = 10;
egm = data(idx).signal{1,1};

samplingFreq = 2*500;
egm = resample(egm,2,1);
plotSignalAndSpectrum(egm, samplingFreq);


% % Use FFT
% Y = fft(egm);
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f = samplingFreq*(0:(L/2))/L;
% figureRight,
% plot(f,P1)
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% 
% 
% 
% figureLeft,
% plot(egm)


[powSpectDens, psdOmega] = ...
  getPowerSpectralDensityBotteron(egm, samplingFreq);

figureRight,
plot(psdOmega, powSpectDens)

%%

plotAnalysisOverMeshInteractive(meshData, cartoData,cartoPointAnalysis, botteronFrequency );



%%

absEgm = abs(egm);
figureLeft,
plot(absEgm)
egmStats = Utils_basic_stats(absEgm);

RPeakIdx = find(absEgm > (egmStats.Mu + egmStats.Std));
figureRight,
plot(egm, 'r')
hold on
plot(RPeakIdx, egm(RPeakIdx), 'o')




%%
% figureRight,
% plot(pointDescriptor);




%% Compare to the full mesh

figure,
trisurf(facets ,meshVertex(:,1), meshVertex(:,2), meshVertex(:,3),...
  egmData(timeFrame, :));


%% Compare to the full mesh

% Save screenshots of the mesh across time

%[currentAL, currentEL] = view(gca);

currentAL = -9.1;
currentEL = -17.2;
timeFrame = 10:200:5002;

resultsFolder = '/Users/carlosAguilar/Documents/Matlab PhD/git synced/cartoAnalysis/simulation';

figResultsFolder = fullfile(resultsFolder);
if ~isdir(figResultsFolder)
    mkdir(figResultsFolder);
end

numSlices = numel(timeFrame);

figure,

for iSlice = 1:numSlices

trisurf(facets ,meshVertex(:,1), meshVertex(:,2), meshVertex(:,3),...
  egmData(timeFrame(iSlice), :));
figTitle = ['TimeSlice ', num2str(timeFrame(iSlice))];
title(figTitle);
figPathName = fullfile(figResultsFolder, [figTitle, '.png']);
hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
%close(currentFig);
fprintf('\n');

end

%%
% figure,
% plot(pointDescriptor, 'LineWidth', 2.5);
% hold on
% plot(egmData(timeFrame, :), 'r');


%cartoAnalysis = load(currentMatFile)
%plotAnalysisOverMeshInteractive(meshData, cartoData,cartoPointAnalysis, botteronFrequency );
