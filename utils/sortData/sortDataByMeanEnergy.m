function [inputDataSorted, meanEnergyValues, idxNRG] = sortDataByMeanEnergy(inputData)

% SORTDATABYMEANENERGY sort data in descending mean energy order.
%

  meanEnergy = nanmean(inputData.^2, 2);
  [meanEnergyValues, idxNRG] = sort(meanEnergy, 'descend');
  inputDataSorted = inputData(idxNRG, :)  ;    
    
% TO-DO: Add here code for svn/gitHub  
% $Revision: $
% $Author: $
% $Date: $
% Copyright