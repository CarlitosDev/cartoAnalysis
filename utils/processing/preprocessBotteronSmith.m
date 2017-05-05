function botteronPreprocessed = preprocessBotteronSmith(inputSignal, samplingFreq)

%% PREPROCESSBOTTEROM preprocess ECG data as described by Botteron-Smith.
%
% Botteron-Smith: A technique for measurement of the extent of spatial
% organization of atrial activation during atrial fibrillation in the
% intact human heart.
% 
% http://www.ncbi.nlm.nih.gov/pubmed/7790014
%
% 
%
% dataIn must be a struct with the fields:
%   dataIn.signalNames
%   dataIn.signalData
%   samplingFreq
% 
% Notes: 
%   In the paper they use a "digital, zero-phase, third-order Butterworth
%   filter". We'll use FIR instead and remove the mean after the
%   processing.
%
% See also filterChebyShev.m
%
% Carlos Aguilar - November 14th 2k15
% Last update Feb 2k16

% TO-DO: Get this paper
% http://www.pubfacts.com/detail/24219142/On-the-preprocessing-of-atrial-electrograms-in-atrial-fibrillation-understanding-Botterons-approach


% Band-passed filtered using a digital, zero-phase, third-order Butterworth
% filter with cutoffs of 40-250 Hz. The absolute value of the output of the
% bandpass filter was then low-pass filtered using a similar third-order
% Butterworth filter with a 20-Hz cut-off

filterType = 'butterworth';
%filterType = 'chebyshev';
%filterType = 'fir';

f1HighPassCutOff = 40 ;
f1LowPassCutOff  = 250;
f2LowPassCutOff  = 20 ;


switch(lower(filterType))
    
    case 'butterworth'

      filterOrder = 3;

      [filtNum, filtDenom] = butter(filterOrder, ...
                              f1HighPassCutOff/samplingFreq, 'high');
      dataSignalProc       = filtfilt(filtNum, filtDenom, inputSignal);

      [filtNum, filtDenom] = butter(filterOrder, ...
                              f1LowPassCutOff/samplingFreq, 'low');
      dataSignalProc       = filtfilt(filtNum, filtDenom, dataSignalProc);

      dataSignalProc       = arrayfun(@abs, dataSignalProc);

      [filtNum, filtDenom] = butter(filterOrder, ...
                              f2LowPassCutOff/samplingFreq, 'low');
      dataSignalProc       = filtfilt(filtNum, filtDenom, dataSignalProc);                           

                            
    case 'chebyshev'
        
     dataSignalProc = filterChebyShev(inputSignal, samplingFreq, ...
        'filterOrder'     , 6    , ...
        'filterRipple'    , 40   , ...
        'filterFrequency' , f1HighPassCutOff , ...
        'filterType'      ,'high', ...
        'showFilterDesign', false);   %#ok<*UNRCH>
    
    dataSignalProc = filterChebyShev(dataSignalProc, samplingFreq, ...
        'filterOrder'     , 6    , ...
        'filterRipple'    , 40   , ...
        'filterFrequency' , f1LowPassCutOff , ...
        'filterType'      ,'low' , ...
        'showFilterDesign', false);  
    
	 dataSignalProc = arrayfun(@abs, dataSignalProc);
    
     dataSignalProc = filterChebyShev(dataSignalProc, samplingFreq, ...
        'filterOrder'     , 6    , ...
        'filterRipple'    , 40   , ...
        'filterFrequency' , f2LowPassCutOff , ...
        'filterType'      ,'low' , ...
        'showFilterDesign', false);         
      
    case 'fir'
        
        dataSignalProc = filterFIR(inputSignal, samplingFreq, ...
                'filterOrder'     , 256, ...
                'filterType'      ,'high', ...
                'filterFrequency' , f1HighPassCutOff, ...
                'showFilterDesign', false);

        dataSignalProc = filterFIR(dataSignalProc, samplingFreq, ...
                'filterOrder'     , 256, ...
                'filterType'      ,'low', ...
                'filterFrequency' , f1LowPassCutOff, ...        
                'showFilterDesign', false);   

        dataSignalProc = arrayfun(@abs, dataSignalProc);

        dataSignalProc = filterFIR(dataSignalProc, samplingFreq, ...
                'filterOrder'     , 256, ...
                'filterType'      ,'low', ...
                'filterFrequency' , f2LowPassCutOff, ...        
                'showFilterDesign', false);  
end


    
% remove the mean
dataSignalProcWhite = ...
    bsxfun(@minus, dataSignalProc, mean(dataSignalProc,2));

botteronPreprocessed = dataSignalProcWhite;

% figure,
% plot(inputSignal, 'LineWidth', 1.5, 'Color', [1 0 0]);
% hold on
% plot(dataSignalProcWhite, 'LineWidth', 1.5, 'Color', [0 1 0]);
