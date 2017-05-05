function dataFromExportFiles = readECGExportFolder(dataFolder)


filePreamble = '2-AI_P';
chunkSize    = 11;

%% scan the folder
dataFileList = dir(fullfile(dataFolder, '*ECG_Export.txt'));


% scan the first file to get the import options
idxFile = 1;
currentFile = fullfile(dataFolder, dataFileList(idxFile).name);

% autoscan file contents (I guess all files look the same)
opts = detectImportOptions(currentFile, ...
  'NumHeaderLines', 3, ...
  'FileType', 'text');

numFiles            = numel(dataFileList);
% compare the new column with the previous one to decide whether importing
% or nor
prevFirstColumn = [];
idxInternal     = 0;
info = [];

for idxFile = 1:2%numFiles
  
  fprintf('Processing %d/%d\n', idxFile, numFiles);
  
  currentFile = fullfile(dataFolder, dataFileList(idxFile).name);
  %opts.VariableNames
  
  % read the table and compare the first column to the one from the
  % previous file
  T = readtable(currentFile, opts);
  if idxFile == 1
    prevFirstColumn = T{:, 1};
    columnsDiff = 1;
  else
    columnsDiff = sum(abs(prevFirstColumn-T{:, 1}));
  end
  
  if columnsDiff > 0
    
    fileID = fopen(currentFile);
    fileHeaders = textscan(fileID, '%s', 4, 'Delimiter', '\n');
    fclose(fileID);
        
    currentID  = extractBetween(dataFileList(idxFile).name, filePreamble, '_');
    
    gainInfo   = fileHeaders{1,1}{2};
    signalGain = str2double(extractAfter(gainInfo,'= '));
    
    mappingInfo = fileHeaders{1,1}{3};
    refUnipolar = extractBetween(mappingInfo, 'Unipolar Mapping Channel=', ...
                  ' Bipolar Mapping Channel=');
    refBipolar  = extractBetween(mappingInfo, 'Bipolar Mapping Channel=', ...
                  ' Reference Channel');
    refChannel  = extractAfter(mappingInfo, 'Reference Channel=');
    
    signalsInfo = fileHeaders{1,1}{4};
    signalNames = split(signalsInfo, ' ');
    signalNames(signalNames == '') = [];        

    idxInternal = idxInternal + 1;
    info(idxInternal).ID           = str2double(currentID);    
    info(idxInternal).gain         = signalGain;
    info(idxInternal).ref_Unipolar = char(refUnipolar);
    info(idxInternal).ref_Bipolar  = char(refBipolar);
    info(idxInternal).ref_Channel  = char(refChannel);
    info(idxInternal).tipo_ECG     = cellstr(signalNames);  
    info(idxInternal).data         = table2cell(T);    
  else
    fprintf('Discarding file %s\n', dataFileList(idxFile).name);
  end
  
  if idxInternal == chunkSize
    fprintf('Saving mat file...');
    save(fullfile(dataFolder, ['ECG', getTimeStamp, '.mat']), 'info');
    info = [];
    idxInternal = 0;
    fprintf('done.\n')
  end
  
end

dataFromExportFiles = info;