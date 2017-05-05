function [inputDataSorted, eKurtValues, idxKurt] = sortDataByKurtosis(inputData)

% SORTDATABYKURTOSIS sort data in ascending excess kurtosis order.
%
%

    % fit the control point to a cubic spline
    % work this out...
    statsTbxLicensed = license('test', 'Statistics_Toolbox');
    if statsTbxLicensed    
        eKurtosis = kurtosis(inputData, 1, 2) - 3;    
    else
        % Need to tile the output of nanmean to center X.
        dim = 2;
        tile = ones(1,max(ndims(inputData),dim));
        tile(dim) = size(inputData,dim);

        % Center X, compute its fourth and second moments, and compute the
        % uncorrected kurtosis.
        x0 = inputData - repmat(nanmean(inputData,dim), tile);
        s2 = nanmean(x0.^2,dim); % this is the biased variance estimator
        m4 = nanmean(x0.^4,dim);
        eKurtosis = (m4 ./ s2.^2) - 3;
    end


    [eKurtValues, idxKurt] = sort(eKurtosis) ;
    inputDataSorted = inputData(idxKurt, :)  ;    
    
% TO-DO: Add here code for svn/gitHub  
% $Revision: $
% $Author: $
% $Date: $
% Copyright