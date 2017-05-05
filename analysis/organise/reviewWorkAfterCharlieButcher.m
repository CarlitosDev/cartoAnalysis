
inputSignal = selDataIn;
samplingFreq = 1e3

figure,
plot(inputSignal)
hold on
plot(dataForICA, 'r')



plotSignalAndSpectrum(inputSignal, samplingFreq)
plotSignalAndSpectrumBotteron(inputSignal, samplingFreq)

%% FIR filtering

samplingFreq = 1e3;
  
    
dataSignalProc = filterFIR(inputSignal, samplingFreq, ...
                'filterOrder'     , 512, ...
                'filterType'      ,'high', ...
                'filterFrequency' , 1.0, ...
                'showFilterDesign', false);
        
figure,
plot(inputSignal)
hold on
plot(dataSignalProc, 'r')

% Apply a notch filter

d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
               'DesignMethod','butter','SampleRate',samplingFreq);


dataSignalProc = filtfilt(d, dataSignalProc);


figure,
plot(inputSignal)
hold on
plot(dataSignalProc, 'r')


dataSignalProc = filterFIR(dataSignalProc, samplingFreq, ...
                'filterOrder'     , 512, ...
                'filterType'      ,'low', ...
                'filterFrequency' , 250.0, ...
                'showFilterDesign', false);
figure,
plot(inputSignal)
hold on
plot(dataSignalProc, 'r') 

%% FIR filtering

samplingFreq = 1e3;
  
    
dataSignalProc = filterFIR(inputSignal, samplingFreq, ...
                'filterOrder'     , 512, ...
                'filterType'      ,'high', ...
                'filterFrequency' , 1.0, ...
                'showFilterDesign', false);
        
figure,
plot(inputSignal)
hold on
plot(dataSignalProc, 'r')

% Apply a notch filter

d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',49,'HalfPowerFrequency2',51, ...
               'DesignMethod','butter','SampleRate',samplingFreq);


dataSignalProc = filtfilt(d, dataSignalProc);


figure,
plot(inputSignal)
hold on
plot(dataSignalProc, 'r')


dataSignalProc = filterFIR(dataSignalProc, samplingFreq, ...
                'filterOrder'     , 512, ...
                'filterType'      ,'low', ...
                'filterFrequency' , 250.0, ...
                'showFilterDesign', false);
figure,
plot(inputSignal)
hold on
plot(dataSignalProc, 'r') 

%% IIR filtering Schilling (page 145)

dataSignalProc = filterFIR(selDataIn, samplingFreq, ...
                'filterOrder'     , 512, ...
                'filterType'      ,'high', ...
                'filterFrequency' , 30.0, ...
                'showFilterDesign', false);


dataSignalProc = filterFIR(dataSignalProc, samplingFreq, ...
                'filterOrder'     , 512, ...
                'filterType'      ,'low', ...
                'filterFrequency' , 150.0, ...
                'showFilterDesign', false);


plotSignalAndSpectrum(dataSignalProc, dataIn.fs )


figure,
plot(inputSignal)
hold on
plot(dataSignalProc, 'r')

%%

%path2CartoFile = fullfile('data','Auriculas','#4', 'Carto1', 'EG.mat');
path2CartoFile = fullfile('data','Auriculas','#2', 'Carto', 'EG.mat');
fullPath2CartoFile = fullfile(rootFolder(), path2CartoFile); 
    
listOfPoints = 116;
lisftOfSignals = {'M1_MINUS_M2'};

        
    variables.fs = 1000;   
    
    % load CartoXP data
[data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile);
  
idx      = find(pointsInfo == listOfPoints);
ecgData  = data(idx);
    

fprintf('\nProcessing point %d\n', listOfPoints);

% Preprocess input data. Be verbose
dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;

% for idx=1:numel(dataIn.signalNames)
%     fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
% end

dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;

[signalPresent, idxSignals] = ismember(lisftOfSignals, dataIn.signalNames);
assert(all(signalPresent), 'Missing signals');
selDataIn = dataIn.signalData(idxSignals, :);

plotSignalAndSpectrum(selDataIn, dataIn.fs )

%%
botteronOut = preprocessBotteronSmith (dataIn.signalData, dataIn.fs );
botteronSignal = botteronOut(idxSignals, :);

plotSignalAndSpectrum(botteronSignal, dataIn.fs )


%%

% Apply Bottteron-Smith using FIR filters


    dataSignalProc = filterFIR(selDataIn, dataIn.fs, ...
            'filterOrder'     , 256, ...
            'filterType'      ,'high', ...
            'filterFrequency' , f1HighPassCutOff, ...
            'showFilterDesign', false);
    
        plotSignalAndSpectrum(dataSignalProc, dataIn.fs )

    dataSignalProc = filterFIR(dataSignalProc, dataIn.fs, ...
            'filterOrder'     , 256, ...
            'filterType'      ,'low', ...
            'filterFrequency' , f1LowPassCutOff, ...        
            'showFilterDesign', false);   
        
           plotSignalAndSpectrum(dataSignalProc, dataIn.fs )
           
    dataSignalProc = arrayfun(@abs, dataSignalProc);
        
      plotSignalAndSpectrum(dataSignalProc, dataIn.fs )
    
    dataSignalProc = filterFIR(dataSignalProc, dataIn.fs, ...
            'filterOrder'     , 256, ...
            'filterType'      ,'low', ...
            'filterFrequency' , f2LowPassCutOff, ...        
            'showFilterDesign', false);
        
     plotSignalAndSpectrum(dataSignalProc, dataIn.fs )  
     
dataSignalProc = bsxfun(@minus, dataSignalProc, mean(dataSignalProc,2));
        
             plotSignalAndSpectrum(dataSignalProc, dataIn.fs ) 
     
%%
        
% Apply Bottteron-Smith using IIR filters

dataSignalProc = filterChebyShev(selDataIn, dataIn.fs, ...
    'filterOrder'     , 6    , ...
    'filterRipple'    , 25   , ...
    'filterFrequency' , f1HighPassCutOff , ...
    'filterType'      ,'high', ...
    'showFilterDesign', false);  

       plotSignalAndSpectrum(dataSignalProc, dataIn.fs )