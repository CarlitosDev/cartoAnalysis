
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

samplingFreq = 1e3;
inputEGM     = dataForICA;
inputECG     = referenceSignal;

plotEGMandECGTimeAndSpectrum(inputEGM, inputECG, samplingFreq);



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