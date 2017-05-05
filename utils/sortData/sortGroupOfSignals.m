function [sortingValues, idxSorted] = sortGroupOfSignals(dataIn, signalsToResearch, varargin)
% SORTGROUPOFSIGNALS using a criteria based on mean energy, kurtosis or PSD
% bands
%
%
% Carlos Aguilar - November 11th 2k16


%% Process input arguments.
    p = inputParser;
    p.FunctionName = 'sortGroupOfSignals';

    vStruct       = @(x) isa(x, 'struct');
    vChar         = @(x) isa(x, 'char');
    vCellString   = @(x) iscellstr(x);
    vSignals2Look = @(x) vChar(x) || vCellString(x);
    vCriteriaType = @(x) find(strcmp(x, {'meanEnergy', 'kurtosis', 'psdBands'}));

    p.addRequired  ('dataIn'           ,               vStruct);
    p.addRequired  ('signalsToResearch',               vSignals2Look);
    p.addOptional  ('sortingCriteria'  , 'meanEnergy', vCriteriaType);


    p.parse(dataIn, signalsToResearch, varargin{:});

    dataIn            = p.Results.dataIn;
    signalsToResearch = p.Results.signalsToResearch;
    sortingCriteria   = p.Results.sortingCriteria;
  
%% Sort the data
 
 [signalPresent, idxSignals] = ismember(signalsToResearch, dataIn.signalNames);
 assert(all(signalPresent), 'Missing signals');
 
 inputData = dataIn.signalData(idxSignals, :);
 
 switch(sortingCriteria)
     case 'meanEnergy'
         [~, sortingValues, idxSorted] = sortDataByMeanEnergy(inputData);
     case 'kurtosis'
         [~, sortingValues, idxSorted]  = sortDataByKurtosis(inputData);
     case 'psdBands'
         warning('Bands are hard coded to (10,30) Hz');
         lowerFreq = 10;
         upperFreq = 30;
         [~, sortingValues, idxSorted] = ...
             sortDataByPSDbands(inputData, samplingFrequency, lowerFreq, upperFreq);
 end
 
 idxSorted = idxSignals(idxSorted);