function atrialFibrillation = synthesiseAtrialFibrillation(numSamples, samplingFreq, varargin)

%% 
% 
% Synthesize an Atrial Fibrillation signal as described in ~\cite{AFPreproc:stridh2001AFCharact, AFPreproc:alcaraz2009sample}
% 
%   main equations. TO-DO: LaTeX them
%
%   theta = (2*pi*n*F0/Fs) + (deltaF/Ff)*sin(n*Ff/fs);
%   ak    = (2/(pi*k)) *(a + deltaA*sin(2*pi*Fa*n/Fs));
%   af(n) = -sum( ak*sin(k*theta));
% 
% Input parameters:
% 
%   numSamples:       Number of samples for the AF signal. 
%   samplingFreq:     Sampling frequency.
%   numHarmonics:     Number of harmonics. Default 15.
%   mainFrequency:    Fundamental frequency of the fibrillation waveform. Default 6 Hz.
%   maxFreqDeviation: Maximum frequency deviation from mainFrequency. Default 3 Hz.  
%   modFrequency:     Modulation frequency of the fibrillatory wave. Default 4 Hz.
%   afAmplitude:      Amplitude. Default 10 uV.
%   modPeakAmplitude: Modulation peak amplitude. Default 10 uV;
%   modFreqAmplitude: Amplitude for the modulation frequency. Default 9 Hz.
% 
% ~\cite{AFPreproc:stridh2001AFCharact}
% Stridh et al.
% Characterization of atrial fibrillation using the surface ECG: time-dependent spectral properties
% IEEE Transactions on Biomedical Engineering 2001
% 
% ~\cite{AFPreproc:alcaraz2009sample}
% Alcaraz et al.
% Medical engineering & physics. 2009
% Sample entropy of the main atrial wave predicts spontaneous termination of paroxysmal atrial fibrillation
%
% 
% Usage:
%   ....
%
% Carlos Aguilar - June 25th 2k15


% Process input arguments.
  p = inputParser;
  p.FunctionName = 'synthesizeAtrialFibrillation';

  vDouble     = @(x) isa(x, 'double');

  p.addRequired  ('numSamples'      ,      vDouble);
  p.addRequired  ('samplingFreq'    ,      vDouble);
  p.addOptional  ('numHarmonics'    ,  15, vDouble);
  p.addOptional  ('mainFrequency'   ,  6 , vDouble);
  p.addOptional  ('maxFreqDeviation',  3 , vDouble);
  p.addOptional  ('modFrequency'    ,  4 , vDouble);
  p.addOptional  ('afAmplitude'     ,  10, vDouble);
  p.addOptional  ('modPeakAmplitude',  10, vDouble);
  p.addOptional  ('modFreqAmplitude',  9 , vDouble);
   
  p.parse(numSamples, samplingFreq, varargin{:});

  
  % Following notation from ~\cite{AFPreproc:alcaraz2009sample}  
    numSamples = p.Results.numSamples;
    M          = p.Results.numHarmonics;
    Fs         = p.Results.samplingFreq;
    F0         = p.Results.mainFrequency;
    deltaF     = p.Results.maxFreqDeviation;
    Ff         = p.Results.modFrequency;
    deltaA     = p.Results.modPeakAmplitude;
    a          = p.Results.afAmplitude;
    Fa         = p.Results.modFreqAmplitude;

  
  
    atrialFibrillation = zeros(1, numSamples);

    for n=1:numSamples
        theta = (2*pi*n*F0/Fs) + (deltaF/Ff)*sin(n*Ff/Fs);    
        for k=1:M
            ak    = (2/(pi*k)) *(a + deltaA*sin(2*pi*Fa*n/Fs));        
            atrialFibrillation(n) = atrialFibrillation(n) - (ak*sin(k*theta));
        end
    end
