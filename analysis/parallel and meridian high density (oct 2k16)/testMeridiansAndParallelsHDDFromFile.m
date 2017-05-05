%% testMeridiansAndParallelsHDDFromFile
%
% This script loads the new cartoXP data stored in .\preprocessed
% data\Results FA_RamonYCajal ECG.mat
% 
%

%%

%currentTestFolder = rootFolder();

preprocessedData = load('D:\carlosAguilar\PhD Carlos\Data (ungitted)\preprocessed data\Results FA_RamonYCajal ECG.mat.mat');


%%
samplingFreq = 1e3;
currentIdx = 1;
inputEGM   =  preprocessedData.cartoPointAnalysis{currentIdx}.dataForICA;
inputECG   =  preprocessedData.cartoPointAnalysis{currentIdx}.referenceSignal;
plotEGMandECGTimeAndSpectrum(inputEGM, inputECG, samplingFreq);

plotCartoAnalysisOverlay(preprocessedData.cartoPointAnalysis{currentIdx});

%% Check the new visualisation
%plotMeshAndAllEGMAnalysis( meshData, cartoData, cartoPointAnalysis, botteronFrequency)
plotMeshAndAllEGMAnalysisParallel( preprocessedData.meshData, preprocessedData.cartoData, ...
  preprocessedData.cartoPointAnalysis, preprocessedData.botteronFrequency)