

lisftOfSignals = 'R2';


[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);

inputSignal = dataIn.signalData(idxSignals, :);

samplingFreq = dataIn.fs;

plotSignalAndSpectrum(inputSignal, samplingFreq);


% filter

filterOrder    = 128;
dataSignalProc = filterFIR(inputSignal, dataIn.fs, ...
        'filterOrder'     , filterOrder, ...
        'filterType'      ,'high'      , ...
        'filterFrequency' , 3.0       , ...
        'showFilterDesign', false);
    
dataSignalProc = filterFIR(dataSignalProc, dataIn.fs, ...
        'filterOrder'     , filterOrder, ...
        'filterType'      ,'low'       , ...
        'filterFrequency' , 30         , ...        
        'showFilterDesign', false);

plotSignalAndSpectrum(dataSignalProc, samplingFreq);