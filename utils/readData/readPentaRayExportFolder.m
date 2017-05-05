function dataFromExportFiles = readPentaRayExportFolder(dataFolder, varNames, ...
  filePreamble, electrodeFileName)

% This function reads the data exported from PentaRay

  chunkSize    = 2500;

  if strcmpi(varNames, 'all')
    allVarnames = true;
  else
    allVarnames = false;
  end

%% Scan the folder and get the options to import the ECG files

  dataFileList = dir(fullfile(dataFolder, '*ECG_Export.txt'));

  % get filename numbers
  fileNums = str2double(extractBetween({dataFileList.name}, '_P', '_'));
  [sortedFileNums, sortedIdx] = sort(fileNums);

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


%% Also import the electrodes on annotation file

  electrodesFileList = dir(fullfile(dataFolder, ['*', electrodeFileName]));
  electrodesFile     = electrodesFileList(idxFile).name;
  electrodesID       = extractBetween(electrodesFile, filePreamble(1:end-1), '_');
  electrodesFilePath = fullfile(dataFolder, electrodesFile);

  % autoscan file contents (I guess all files look the same)
  electrodesOpts = detectImportOptions(electrodesFilePath, ...
    'NumHeaderLines', 1, ...
    'FileType', 'text');


%%


for idxFile = 1:numFiles
  
  internalIdx     = sortedIdx(idxFile);
  currentFilename = dataFileList(internalIdx).name;
  currentFile     = fullfile(dataFolder, currentFilename);
  
	fprintf('Processing %s (%d/%d)\n', currentFilename, idxFile, numFiles);
  
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
    
    prevFirstColumn = T{:, 1};
    
    fileID      = fopen(currentFile);
    fileHeaders = textscan(fileID, '%s', 4, 'Delimiter', '\n');
    fclose(fileID);
        
    currentID  = extractBetween(currentFilename, filePreamble, '_');
    
    gainInfo   = fileHeaders{1,1}{2};
    signalGain = str2double(extractAfter(gainInfo,'= '));
    
    mappingInfo = fileHeaders{1,1}{3};
    refUnipolar = extractBetween(mappingInfo, 'Unipolar Mapping Channel=', ...
                  ' Bipolar Mapping Channel=');
    refBipolar  = extractBetween(mappingInfo, 'Bipolar Mapping Channel=', ...
                  ' Reference Channel');
    refChannel  = extractAfter(mappingInfo, 'Reference Channel=');
 
    tableVarNames = T.Properties.VariableNames;
    if allVarnames
      idxVarsToKeep      = true(1, numel(tableVarNames));
      idxVarsToKeep(end) = false;
    else
      idxVarsToKeep = contains(tableVarNames, varNames);
    end
    
    % read the electrodes position
    currentElectrodeFile = ...
      strrep(electrodesFile, electrodesID, ['P', num2str(sortedFileNums(idxFile))]);
    electrodesFilePath   = fullfile(dataFolder, currentElectrodeFile);
    electrodesTable      = readMagneticElectrodes(electrodesFilePath, electrodesOpts);
    
    idxInternal = idxInternal + 1;
    info(idxInternal).ID            = str2double(currentID);
    info(idxInternal).gain          = signalGain;
    info(idxInternal).ref_Unipolar  = char(refUnipolar);
    info(idxInternal).ref_Bipolar   = char(refBipolar);
    info(idxInternal).ref_Channel   = char(refChannel);
    info(idxInternal).tipo_ECG      = cellstr(tableVarNames(idxVarsToKeep));
    % to match Marga's format...
    info(idxInternal).data          = num2cell(T{:, idxVarsToKeep}', 2)';
    info(idxInternal).electrodesPos = electrodesTable;
  else
    fprintf('Discarding file %s\n', currentFilename);
  end
  
  if idxInternal == chunkSize
    fprintf('Saving mat file...');
    save(fullfile(dataFolder, ['ECG', getTimeStamp(), '.mat']), 'info');
    info = [];
    idxInternal = 0;
    fprintf('done.\n')
  end
  
end

dataFromExportFiles = info;