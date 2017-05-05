function [powSpectDens, psdOmega] = getPowerSpectralDensityBotteron(inputSignal, samplingFreq, varargin)


%% GETPOWERSPECTRALDENSITYBOTTERON get PSD for atrial EGMs based on Botteron-Smith's paper.
%
% Preprocess data based on Botteron-Smith's paper for EGM frequency
% analysis and get PSD using Welch-WOSA modified periodogram.
%
% 
% Usage
%   filteredData = getPowerSpectralDensityBotteron(inputSignal, samplingFreq, ...
%    'windowSize'   , 4096, ...
%    'lengthFFT'    , 8192, ...
%    'psdCutOff'    , 0.9 , ...
%    'numberOverlap', 2048);
%
% Carlos Aguilar - November 14th 2k15



% Process input arguments.
  p = inputParser;
  p.FunctionName = 'getPowerSpectralDensityBotteron';

  vDouble     = @(x) isa(x, 'double');
  vWindow     = @(x) isa(x, 'char') & any(strcmpi(x, {'hamming', 'rectangular'}));
  getNFFTSz   = @(x) 2^nextpow2(length(x));

  p.addRequired('inputSignal'  ,                          vDouble);
  p.addRequired('samplingFreq' ,                          vDouble);
  p.addOptional('windowType'   ,            'hamming'  ,  vWindow);
  p.addOptional('windowSize'   ,                 4096  ,  vDouble);
  p.addOptional('lengthFFT'    , getNFFTSz(inputSignal),  vDouble);
  p.addOptional('numberOverlap',                   []  ,  vDouble);
  p.addOptional('psdCutOff'    ,                 1.00  ,  vDouble);

  p.parse(inputSignal, samplingFreq, varargin{:});

  inputSignal    = p.Results.inputSignal  ;
  samplingFreq   = p.Results.samplingFreq ;
  windowType     = p.Results.windowType   ;
  windowSize     = p.Results.windowSize   ;
  lengthFFT      = p.Results.lengthFFT    ;
  numberOverlap  = p.Results.numberOverlap;
  psdCutOff      = p.Results.psdCutOff    ;

  % Apply Botteron-Smith preprocessing to the dataset.
  inputSignal = preprocessBotteronSmith(inputSignal, samplingFreq);
  
  
  % Call the built-in method to get the PSD.
  if strcmpi(windowType, 'hamming')
  [pxx, wxx] = ...
    pwelch(inputSignal, windowSize, numberOverlap, lengthFFT, ...
    samplingFreq, 'onesided');
  else
      inputSignalLenght = length(inputSignal);
  [pxx, wxx] = periodogram(inputSignal, rectwin(inputSignalLenght), ...
    inputSignalLenght, samplingFreq);
  end
  
  % cut-off PSD to get LT PSDCUTOFF of the energy
  if psdCutOff < 1.0

    normPxx    = pxx / sum(pxx) ;
    cumNormPxx = cumsum(normPxx);
    validIdx   = cumNormPxx < psdCutOff;
   
    pxx = pxx(validIdx);
    wxx = wxx(validIdx);
    
  elseif  (psdCutOff > 1.0) &&  (psdCutOff <= samplingFreq/2 )

    validIdx = wxx <= psdCutOff;
    pxx      = pxx(validIdx);
    wxx      = wxx(validIdx);
  end
    
  
% output  
powSpectDens = pxx;
psdOmega     = wxx;%wxx*samplingFreq/max(wxx);

    
%    figure, plot(omegaCapped, pxxCapped)    
  
% TO-DO: Add here code for svn/gitHub  
% $Revision: $
% $Author: carlosAguilar $
% $Date: $
% Copyright 2015 OGTel LTD.