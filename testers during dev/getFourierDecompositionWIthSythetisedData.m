%% Generate a 10Hz and 25Hz signal with noise and baseline drift
n   = 1:2500;
Fs  = 1e3;
samplingFreq = Fs;

% (A)
atrialFibrillation = synthesizeAtrialFibrillation(numel(n), samplingFreq, ...
  'mainFrequency', 14);
x1 = atrialFibrillation;

% (B)
f2 = 28;
omega = 2*pi*f2/Fs;
x2 = 2.5*sin(omega * n);

% (C)
noise = 0.75*rand(1, length(n));

% (D)
omega = 2*pi*1.2/Fs;
baselineWander = 1.0*sin(omega * n);

% power line interference 50Hz
omega = 2*pi*50/Fs;
powerLineInterference = 0.25*sin(omega * n);

synSignal = x1+x2+noise+baselineWander+powerLineInterference;

figureRight,
plot(n, synSignal)
title('Input signal');
originalAxes = gca();

inputSignal  = synSignal;
lengthSignal = numel(inputSignal);


%% let's say only high amplitude waves make it to the skin

%hist(atrialFibrillation, 30)
atrialFibrillation2 = zeros(size(atrialFibrillation));
atrialFibrillationAbs = abs(atrialFibrillation);

validIdx = atrialFibrillationAbs>22.0;
atrialFibrillation2(validIdx) = atrialFibrillation(validIdx)*0.85;


validIdx = atrialFibrillationAbs>9.0 & atrialFibrillationAbs<12.0;
atrialFibrillation2(validIdx) = atrialFibrillation(validIdx)*0.40;

validIdx = atrialFibrillationAbs>6.0 & atrialFibrillationAbs<9.0;
atrialFibrillation2(validIdx) = atrialFibrillation(validIdx)*0.25;

validIdx = atrialFibrillationAbs>4.0 & atrialFibrillationAbs<6.0;
atrialFibrillation2(validIdx) = atrialFibrillation(validIdx)*0.15;

validIdx = atrialFibrillationAbs>2.0 & atrialFibrillationAbs<4.0;
atrialFibrillation2(validIdx) = atrialFibrillation(validIdx)*0.1;

validIdx = atrialFibrillationAbs>0.0 & atrialFibrillationAbs<2.0;
atrialFibrillation2(validIdx) = atrialFibrillation(validIdx)*0.05;

noise = 0.95*rand(1, length(n));

% power line interference 50Hz
omega = 2*pi*50/Fs;
powerLineInterference = 0.25*sin(omega * n);

synSignal2 = atrialFibrillation2+noise+powerLineInterference;

figureRight,
plot(n, synSignal, 'r', 'LineWidth', 1.0)
hold on
plot(n, synSignal2, 'LineWidth', 1.5);
title('Input signal');


%% get the FS representation

[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(synSignal, samplingFreq, ...
                            'lengthFFT', 512, 'doCoeffsPlot', false);
                          
%data = synSignal;               
%%
data = synSignal'; 
nmbOfCoeffs = 80;
mFileName = 'as'
[FS_data,Ck,a0,an,bn] = FSAproximationNL ( data, nmbOfCoeffs, mFileName );
                          
                          %%

                          
  [yHatS2, fourierSeriesInfoS2] = getSignalAsFourierSeries(synSignal2, samplingFreq, ...
    'lengthFFT', 512, 'doCoeffsPlot', true);


  %%


samplingFreq = 1e3;

inputEGM = dataForICA;
inputECG = referenceSignal;

data = inputEGM;
Fs = samplingFreq;

inputSignal = data;

plotEGMandECGTimeAndSpectrum(inputEGM, inputECG, samplingFreq)


plotSignalAndSpectrum(inputEGM, samplingFreq)

plotSignalAndSpectrumBotteron(inputEGM, samplingFreq)

%%

load('cartoProcessed.mat')

%%



[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq);

[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512);

[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512, 'maxCoeffs', 1);
                          
[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512, 'maxCoeffs', 12);
                          
[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512, 'maxCoeffs', 2, ...
                            'narrowSearch', [2 18]);
                          
[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512, 'maxCoeffs', 3, ...
                            'narrowSearch', [3 10]);             
                          
                          
[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 1024, 'maxCoeffs', 8, ...
                            'narrowSearch', [2 18]);
                          
%%

n_harmonics = 1;
f0 = 7.5;
[y_est,e_est,probs] = myHarmonicSeries(data,Fs,n_harmonics,f0);

plotSignalAndSpectrum(y_est, samplingFreq)

plotEGMandECGTimeAndSpectrum(inputEGM, y_est, samplingFreq)