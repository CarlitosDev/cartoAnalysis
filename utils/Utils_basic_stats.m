function basicStats = Utils_basic_stats( data )

% UTILS_BASIC_STATS get the basic statistics on a dataset


% Flags
plotPDF  = false;
populate = false;

% Get only valid indexes
validIdx     = ~isnan(data);
data         = data(validIdx);
validSamples = numel(data);

basicStats              = [];
basicStats.validSamples = validSamples;


% Create an struct with useful information
%   $\mu$,$\sigma$,median,$Val_min$,$Val_Max$
basicStats.Mu         = mean(data);
basicStats.Std        = std(data);
basicStats.median     = median(data);
basicStats.minVal     = min(data);
basicStats.maxVal     = max(data);
    statsTbxLicensed = license('test', 'Statistics_Toolbox');
    if statsTbxLicensed    
        basicStats.skewness   = skewness(data, 0);
        basicStats.kurtosis   = kurtosis(data, 0);
        basicStats.eKurtosis  = basicStats.kurtosis-3;% excess Kurtosis 
        % Also add caculations for Q1, IQR and Q3
        quantiles             = quantile(data, [0.25, 0.75]);
        basicStats.Q1         = quantiles(1);
        basicStats.Q3         = quantiles(2);
        basicStats.IQR        = quantiles(2) - quantiles(1);
        basicStats.outliersHi = quantiles(2) + 1.5*basicStats.IQR;
        basicStats.outliersLo = quantiles(1) - 1.5*basicStats.IQR;
    else
        basicStats.skewness   = 0;
        basicStats.kurtosis   = 0;
        basicStats.eKurtosis  = 0;% excess Kurtosis
        
        % Also add caculations for Q1, IQR and Q3
        quantiles             = 0;
        basicStats.Q1         = 0;
        basicStats.Q3         = 0;
        basicStats.IQR        = 0;
        basicStats.outliersHi = 0;
        basicStats.outliersLo = 0;
    end

%{
% Look at this...
getBasicStats = @(x) cellfun(@(f) f(x), {@mean, @std, @median, @min, @max});
[basicStats.Mu, basicStats.Std, basicStats.median, basicStats.median, ...
  basicStats.minVal, basicStats.maxVal ]= getBasicStats(currentData);
%}




% Verbose
if populate
    fprintf('\nStats:\n'); 
    fprintf('  Mu        %3.2f\n' , basicStats.Mu );
    fprintf('  Std       %3.2f\n' , basicStats.Std );
    fprintf('  median    %3.2f\n' , basicStats.median );
    fprintf('  minVal    %3.2f\n' , basicStats.minVal );
    fprintf('  maxVal    %3.2f\n' , basicStats.maxVal );
    fprintf('  skewness  %3.2f\n' , basicStats.skewness );
    fprintf('  eKurtosis %3.2f\n' , basicStats.eKurtosis );
    fprintf('  Quant1    %3.2f\n' , basicStats.Q1 ) ;
    fprintf('  Quant3    %3.2f\n' , basicStats.Q3 ) ;
    fprintf('  IQR       %3.2f\n' , basicStats.IQR );
end

% PDF and CDF aproxs
if plotPDF
  figure,
  subplot(121)
  hist(data, 20)
  subplot(122)
  cdfplot( data )
end


% $Revision: $
% $Author: oe68 $
% $Date: $
% Copyright 2015 Tesco PLC.