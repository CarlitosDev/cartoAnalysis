


idxPoint = 1


currentPoint = allPoints(idxPoint);
  
idx      = find(information.npunto == currentPoint);
ecgData  = data(idx);
figTitle =  sprintf('Point %d', currentPoint);

%% Preprocess input data. Be verbose
dataIn = [];
dataIn.signalNames = ecgData.tipo_ECG;

for idx=1:numel(dataIn.signalNames)
    fprintf('%d - %s\n', idx, dataIn.signalNames{idx});
end

dataIn.signalData  = cell2mat(ecgData.signal');
dataIn.fs          = variables.fs;
if doPreprocess
    dataOut = preprocessForICA (dataIn);
else
  dataOut.signalData = dataIn.signalData;
end

%%


i = 1;

fcFilter = [0.5, 55];
dataSignalProc = filterChebyShev(dataIn.signalData, dataIn.fs, ...
        'filterOrder' , 8, 'filterRipple', 30, ...
        'filterFrequency' , fcFilter, ...
        'filterType'      ,'bandpass', ...
        'showFilterDesign', false);
    
    
chebyDesign = sprintf('CHEB N=%d - Fc (%3.2f, %3.2f) Hz', 8, fcFilter(1), fcFilter(2));
    
N = 255+ 1; % por ejemplo
fc = fcFilter;	% Hz
b = fir1(N, fc/dataIn.fs, 'bandpass');	% o con la opción de paso alto
figure,
%freqz(b,1,1024,dataIn.fs);
xfiltrada = filtfilt(b,1,dataIn.signalData(i, :)); % y a correr

firDesign = sprintf('FIR N=%d - Fc (%3.2f, %3.2f) Hz', N, fc(1), fc(2));


figure,
plot(dataIn.signalData(i, :), 'LineWidth',2)
hold on, 
plot(dataSignalProc(i, :),'LineWidth',2, 'Color', 'r')
hold on,
plot(xfiltrada, 'Color', 'g', 'LineWidth',2)

legend('Original', chebyDesign,  firDesign)






%%

i = 5;

fcFilter = [0.1, 55];
dataSignalProc = filterChebyShev(dataIn.signalData, dataIn.fs, ...
        'filterOrder' , 8, 'filterRipple', 30, ...
        'filterFrequency' , fcFilter, ...
        'filterType'      ,'bandpass', ...
        'showFilterDesign', false);
    
    
chebyDesign = sprintf('CHEB N=%d - Fc (%3.2f, %3.2f) Hz', 8, fcFilter(1), fcFilter(2));
    
N = 255+ 1; % por ejemplo
fc = fcFilter;	% Hz
b = fir1(N, fc/dataIn.fs, 'bandpass');	% o con la opción de paso alto
%figure,
%freqz(b,1,1024,dataIn.fs);
xfiltrada = filtfilt(b,1,dataIn.signalData(i, :)); % y a correr

firDesign = sprintf('FIR N=%d - Fc (%3.2f, %3.2f) Hz', N, fc(1), fc(2));


figure,
plot(dataIn.signalData(i, :), 'LineWidth',2)
hold on, 
plot(dataSignalProc(i, :),'LineWidth',2, 'Color', 'r')
hold on,
plot(xfiltrada, 'Color', 'g', 'LineWidth',2)

legend('Original', chebyDesign,  firDesign)

%%



i = 3;


N = 64 + 1; % por ejemplo
fc = 50;	% Hz
b = fir1(N, fc/dataIn.fs);	% o con la opción de paso alto
figure,
freqz(b,1,1024,dataIn.fs);
xfiltrada = filtfilt(b,1,dataIn.signalData(i, :)); % y a correr


figure,
plot(dataIn.signalData(i, :))
hold on, 
plot(dataSignalProc(i, :), 'r')
hold on,
plot(xfiltrada, 'g')