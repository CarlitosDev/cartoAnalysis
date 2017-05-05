
fileName   = 'lpvECG.mat';
fileName   = 'mapaFAECG.mat';
pathToFile = fullfile(getDataFolder(), 'FA_RamonYCajal');
fullFname  = fullfile(pathToFile, fileName);
a   = load(fullFname);
ECG = a.ECG;

%%

fileName   = 'ECG.mat';
pathToFile = fullfile(getDataFolder(), 'RyC_Patient 2016_12_20');
fullFname  = fullfile(pathToFile, fileName);
a   = load(fullFname);
ECG = a.ECG;

%%

% signalNames = unique({ECG(1).tipo_ECG});
signalNames = unique(ECG(1).tipo_ECG);

%%

for idxA = 1:numel(signalNames)
  
  currentSignal = signalNames{idxA};
  idxSymbol     = strfind(currentSignal, '_');
  numChars      = numel(idxSymbol);
  
  if numChars > 3
    signalType = 'BIP';
  else
    signalType = 'MON';
  end
  
  newName = currentSignal(1:idxSymbol(end-1)-1);
  
  fprintf('''%s'', ''%s'', ''%s''; ...\n', currentSignal, ...
    newName, signalType);
  
end
