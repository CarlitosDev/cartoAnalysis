function electrodesTable = readMagneticElectrodes(electrodesFile, electrodesOpts)

electrodesTable = readtable(electrodesFile, electrodesOpts);

varNames = electrodesTable.Properties.VariableNames;
idx = strcmpi('ExtraVar1', varNames);
electrodesTable(:, idx) = [];


electrodeNumbers = array2table([zeros(1,2), 1:height(electrodesTable)-2]', ...
                   'VariableNames',{'electrodeNumbers'});
electrodesTable = [electrodesTable, electrodeNumbers];

