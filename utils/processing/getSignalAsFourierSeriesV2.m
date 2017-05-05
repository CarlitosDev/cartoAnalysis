function [yHat, fourierSeriesInfo] = getSignalAsFourierSeriesV2(inputSignal, samplingFreq, varargin)


%% GETSIGNALASFOURIERSERIES decompose a signal in Fourier's coefficients
%
%
% Carlos Aguilar - May 18th 2k16



  % Process input arguments.
  p = inputParser;
  p.FunctionName = 'getSignalAsFourierSeries';

  vDouble     = @(x) isa(x, 'double');
  vWindow     = @(x) isa(x, 'char') & any(strcmpi(x, {'hamming', 'rectangular'}));
  vLogical    = @islogical;
  getNFFTSz   = @(x) 2^nextpow2(length(x));
  vArraySize2 = @(x) isa(x, 'double') & numel(x) == 2;

  p.addRequired('inputSignal'  ,                          vDouble    );
  p.addRequired('samplingFreq' ,                          vDouble    );
  p.addOptional('lengthFFT'    , getNFFTSz(inputSignal),  vDouble    );
  p.addOptional('windowType'   ,            'hamming'  ,  vWindow    );
  p.addOptional('windowSize'   ,                 4096  ,  vDouble    );
  p.addOptional('findHarmonics',                 true  ,  vLogical   );
  p.addOptional('doTimePlot'   ,                 true  ,  vLogical   );
  p.addOptional('doCoeffsPlot' ,                 false ,  vLogical   );
  p.addOptional('maxCoeffs'    ,                 10    ,  vDouble    );
  p.addOptional('narrowSearch' ,                 []    ,  vArraySize2);

  p.parse(inputSignal, samplingFreq, varargin{:});

  inputSignal    = p.Results.inputSignal  ;
  samplingFreq   = p.Results.samplingFreq ;
  lengthFFT      = p.Results.lengthFFT    ;
  
  % Find 2nd and 3rd harmonics. Needs further checks.
  findHarmonics  = p.Results.findHarmonics;
  
  % flag indicating where to plot
  doTimePlot     = p.Results.doTimePlot;
  doCoeffsPlot   = p.Results.doCoeffsPlot;
  
  % TO-DO: Apply some windowing to smooth the FFT
  windowType     = p.Results.windowType   ;
  windowSize     = p.Results.windowSize   ;  
  
  % Maximum number of coefficients
  maxCoeffs      = p.Results.maxCoeffs    ;  
  
  % Narrow the search for coefficients
  narrowSearch   = p.Results.narrowSearch ;
  
  %% Get the Fourier coefficients through the FFT 
  
  %   FFT
  usedFFT        = lengthFFT/2;
  inputSignalFFT = fft(inputSignal, lengthFFT)/lengthFFT;
  validFFT       = inputSignalFFT(1:usedFFT);
  deltaFFT       = samplingFreq/lengthFFT;
  
  % We can narrow the search by setting a window. To do so, let's set to
  % nil the coeffs outside the window.
  %
  % TO-DO: We can apply some filtering in here to get the Botteron-Smith signal
  % overlap-add might work. 
  %freqAxis      = 0.5*samplingFreq*(1:usedFFT-1)/usedFFT; % ignore 0 Hz (A0)
  freqAxis      = 0.5*samplingFreq*(0:usedFFT)/usedFFT; % ignore 0 Hz (A0)
  %freqAxis      = samplingFreq*(1:usedFFT-1)/usedFFT; % ignore 0 Hz (A0)
  
  if ~isempty(narrowSearch)
    % uncomment the PSD trick if we want to preserve the power
    %currentPDS    = norm(validFFT);
    capFFTIdx     = freqAxis > narrowSearch(1) & freqAxis < narrowSearch(2);
    validFFT(~capFFTIdx) = 0;
    %cappedPDS = norm(validFFT);
    % redistribute the stolen power...
    %validFFT  = validFFT*currentPDS/cappedPDS;
  end

  % Get FS coeffs
  An =  2*real(validFFT);
  Bn = -2*imag(validFFT);
  A0 = An(1)/2;
  An = An(2:end);
  Bn = Bn(2:end);
  
  
  %% Let's get only the prominent coefficients
  
  coefsamplingFreqEnergy  = sqrt(An.^2 + Bn.^2);
  % figure, stem(freqAxis, coefsamplingFreqEnergy)

  % Use FINDPEAKS to do the job. 
  % TO-DO: Use the width of the peaks.
  if license('test', 'Signal_Toolbox')
      [~, idxMainFreqs, mainFreqWidth, mainFreqProm]= ...
          findpeaks(coefsamplingFreqEnergy,'SortStr','descend');
      if numel(idxMainFreqs) > maxCoeffs
          idxMainFreqs = idxMainFreqs(1:maxCoeffs);
      end
  else
      [~, idxMainFreqs] = sort(coefsamplingFreqEnergy, 'descend');
      idxMainFreqs = idxMainFreqs(1:10);
      warning('No Signal_Toolbox. This should be done with method\n');
  end

  mainFreqs = freqAxis(idxMainFreqs);
  numPeaks  = numel(mainFreqs);
  
  if doCoeffsPlot
    figure, 
    %plot([0 freqAxis], [A0.^2 coefsamplingFreqEnergy]);
    plot(freqAxis, coefsamplingFreqEnergy)
    hold on, 
    stem(freqAxis(idxMainFreqs), coefsamplingFreqEnergy(idxMainFreqs));
    title('FFT and prominent Fourier coefficients');
  end


%% Under dev: look for possible 2nd and 3rd harmonics within the peaks
% For every peak that was found, read the remaining ones in the peaks list
% to find a 2nd or 3rd order harmonic. 

  if findHarmonics

    harmonicsIdx = zeros(1, numPeaks);
    harmonicsEps = 0.08;

    for idxPeak = 1:numPeaks-1
        currentPeak       = mainFreqs(idxPeak);
        remainingPeaks    = mainFreqs(1+idxPeak:end);
        possibleHarmonics = remainingPeaks./currentPeak;
        lowerFreqsIdx     = currentPeak > remainingPeaks;
        possibleHarmonics(lowerFreqsIdx) = 1./possibleHarmonics(lowerFreqsIdx);
        nd    = round(possibleHarmonics);
        ndRem = abs(possibleHarmonics-nd);
        harmonics = find(nd == 2 | nd == 3 & ndRem < harmonicsEps);
        % pair the harmonics
        if ~isempty(harmonics) 
            harmonicsIdx(idxPeak)           = idxPeak;
            harmonicsIdx(harmonics+idxPeak) = idxPeak;
        end
    end

    numMatches = max(harmonicsIdx);
    for idx=1:numMatches
        fprintf('\nPair %d of possible harmonics', idx);
        fprintf('\n\t%3.2f Hz', mainFreqs(harmonicsIdx==idx));
    end

  end


%% reconstruct Y hat with the main coefficients

% I'm screwing something up in here...
%%

  An =  2*real(validFFT);
  Bn = -2*imag(validFFT);
  A0 = An(1)/2;
  An = An(2:end);
  Bn = Bn(2:end);

%%

newFFT = zeros(size(validFFT));
newFFT(1) = validFFT(1);



%%
%   FS Coefficients
C0 = real(validFFT(1));
Ck = 2*validFFT(idxMainFreqs);

Ck_Mag = [C0,abs(Ck)];
Ck_Phase = [0.0, angle(Ck)];

numSamples = length(inputSignal);

tR = linspace(0, 1, numSamples);

%yHat  = A0 * ones(1, numSamples);
yHat  = C0 * ones(1, numSamples);

  for k = 1:numPeaks               
    for n = 1:numSamples
          yHat(n) = yHat(n) + real(Ck(k)*exp(1i*2*pi*freqAxis(idxMainFreqs(k))*tR(n)/numSamples));
    end
  end

%%
  numSamples = length(inputSignal);
  
  % $ \omega = fract{2 \pi f}{F_s} $
  %k = (1:usedFFT)*deltaFFT;
  t = 1:numSamples;
  
  %omega = 2*pi*k/samplingFreq;
  omega = 2*pi*freqAxis/samplingFreq;
  
  % Allocate minding the mean
  yHat  = A0 * ones(1, numSamples);

  for k = 1:numPeaks
    coeffIdx   = 1+idxMainFreqs(k);
    yHat       = yHat +  ...
                 An(coeffIdx)*cos(omega(coeffIdx)*t) + ...
                 Bn(coeffIdx)*sin(omega(coeffIdx)*t);
  end

  if doTimePlot
    figure,
    plot(inputSignal, 'LineWidth', 1.0);
    hold on,
    plot(yHat, 'r', 'LineWidth', 1.5);
    fprintf('\nReconstruction with %d coefficients oscillating at\n', numPeaks);
    fprintf('\t%3.2f Hz\n', mainFreqs);
  end
  
  % measures: sum of squared errors of prediction (SSE) and MSE
  sse = sum((yHat-inputSignal).^2);
  mse = sse/numSamples;
  
  %% dump data

  fourierSeriesInfo = [];
  
  fourierSeriesInfo.A0 = A0;
  fourierSeriesInfo.An = An;
  fourierSeriesInfo.Bn = Bn;
  
  fourierSeriesInfo.mainCoeffsIdx = coeffIdx;
  fourierSeriesInfo.omega         = omega;
  fourierSeriesInfo.deltaFFT      = deltaFFT;
  fourierSeriesInfo.mainFreqs     = mainFreqs;
  fourierSeriesInfo.numPeaks      = numPeaks;
  fourierSeriesInfo.mse           = mse;
  fourierSeriesInfo.sse           = sse;
  
  

  % TO-DO: Add here code for svn/gitHub  
  % $Revision: $
  % $Author: carlosAguilar $
  % $Date: $
  % Copyright 2016 OGTel LTD.
  