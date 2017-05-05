%% Generate a 10Hz and 28Hz signal and add some noise and baeline drift
n   = 1:3000;
Fs  = 1e3;

f1  = 10;
f2  = 28;


% (A)
omega = 2*pi*f1/Fs;
x1 = 1.25*sin(omega * n);

% (B)
omega = 2*pi*f2/Fs;
x2 = 1.0*sin(omega * n);

% (C)
noise = 0.65*rand(1, length(n));

% (D)
omega = 2*pi*1.2/Fs;
baselineWander = 1.0*sin(omega * n);

inputSignal = x1+x2+noise+baselineWander;

figure
plot(n, inputSignal)


%%   

getNFFTSz   = @(x) 2^nextpow2(length(x));
lengthFFT   = getNFFTSz(inputSignal)
usedFFT        = lengthFFT/2;
inputSignalFFT = fft(inputSignal, lengthFFT)/lengthFFT;
deltaFFT       = samplingFreq/lengthFFT;
validFFT       = inputSignalFFT(1:usedFFT);

figure,
plot(abs(validFFT))

%%

samplingFreq = Fs;

[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq);

[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 2048);

[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512, 'maxCoeffs', 10);
                          
[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512, 'maxCoeffs', 12);
                          
[yHat, fourierSeriesInfo] = getSignalAsFourierSeries(inputSignal, samplingFreq, ...
                            'lengthFFT', 512, 'maxCoeffs', 2, ...
                            'narrowSearch', [2 18]);
                          