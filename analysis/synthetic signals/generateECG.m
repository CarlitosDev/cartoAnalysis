

matFile = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos\tests\filter benchmark\synthesized data\test1.mat';

figPathName = fullfile(resultsFolder, [figTitle, '.png']);
hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'png');
close(currentFig);
fprintf('\n');



%% Synth ECG + noise
addpath(testFolder(),)
% Call ecgsyn
% IEEE Transactions On Biomedical Engineering, 50(3), 289-294, March 2003.
samplingFreq = 1e3;
numBeats = 12;
Anoise = 0.05;
hrmean = 72;
hrstd = 1;
lfhfratio = 0.5;
sfint = 1000;
ti=[-70 -15 0 15 100];
ai=[1.2 -5 30 -7.5 0.75];
bi=[0.25 0.1 0.1 0.1 0.4];

[simECG, ipeaks, addNoise] = ecgsyn(samplingFreq,numBeats,Anoise,hrmean,hrstd,lfhfratio,sfint,ti,ai,bi);
cleanECG = simECG - addNoise;


%% Synth AF

numSamples = length(simECG);
atrialFibrillation = synthesizeAtrialFibrillation(numSamples, samplingFreq)  ;

% scale AF(uV) to mV
atrialFibrillation = atrialFibrillation/1e3;

% Carlos: Fix the scaling as it doesn't make any sense ?????
atrialFibrillation = 10*atrialFibrillation;

%% Generate ECG with undergoing AF

% add signals
simECGwithAF = atrialFibrillation + simECG;

figure,
subplot(211),
plot(cleanECG),
hold on
plot(addNoise, 'g'),
plot(10*atrialFibrillation, 'r');
legend('synthECG', 'addedNoise', 'synthAF');
title('Synthtesized ECG+AF');
subplot(212),
plot(simECGwithAF)


%% Prepare data for fastICA

% create mixtures
inputSignals = [cleanECG; atrialFibrillation; addNoise];
Amatrix  = rand(size(inputSignals,1));
mixedSig = Amatrix*inputSignals;


figure,
idx = 1;
subplot(3,1,idx), plot(mixedSig(idx, :))
idx = 1 + idx;
subplot(3,1,idx), plot(mixedSig(idx, :), 'r')
idx = 1 + idx;
subplot(3,1,idx), plot(mixedSig(idx, :), 'g')
  

% save data
testData = struct;
testData.simECG             = simECG;
testData.cleanECG           = cleanECG;
testData.addNoise           = addNoise;
testData.atrialFibrillation = atrialFibrillation;
testdata.Amatrix            = Amatrix;
testdata.mixedSig           = mixedSig;

matFile = 'D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\codigoCarlos\tests\filter benchmark\synthesized data\test1.mat';
save(matFile, 'testdata');

