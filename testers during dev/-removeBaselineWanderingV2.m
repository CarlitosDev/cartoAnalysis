function [detrendedData, baselineSignal] = removeBaselineWanderingV2(inputData, samplingPercentage)

% REMOVEBASELINEWANDERING Remove wandering of the base line of an input
% signal.
%
% 
%inputSignal     = inputData;
numberOfSamples = numel(inputData);
extraLength     = round(numberOfSamples*1.2);
extraLength     = 0.5*(extraLength + double(mod(extraLength,2)~=0));
padInputData    = padarray(inputData, [0 extraLength], 'symmetric'); 
padDataLenght   = numel(padInputData);

% 20 percent of the data
stepSize = floor(samplingPercentage*padDataLenght/100);

% correct this shit
% numSmpFactor = factor(numberOfSamples);
% numPrimes    = numel(numSmpFactor);
% idxExcluded  = false(1, numPrimes);
% mostFactors  = numSmpFactor;
% 
% for j=1:numel(numSmpFactor)    
%     idxExcluded(j) = 1;  
%     newFactors = numSmpFactor(j) .* numSmpFactor(~idxExcluded);
%     idxExcluded(j) = 0;
%     mostFactors = [mostFactors, newFactors]; %#ok<AGROW>
% end
% mostFactors = unique(mostFactors);

% numSmpFactor = factor(numberOfSamples);
% [uniqueFactors, ~, ic] = unique(numSmpFactor);
% numberOfChunks = numberOfSamples./uniqueFactors;


%%


% add a snippet to do it automatically.
xControlPoints    = 0:stepSize:padDataLenght;
xControlPoints(1) = 1;

% let's smooth the signal before getting the values for the control points
filterKernel      = ones(1, stepSize)./stepSize;
filteredInputData = conv(padInputData, filterKernel, 'same');
yControlPoints    = filteredInputData(xControlPoints);



% fit the control point to a cubic spline
% work this out...
cfTbxLicensed = license('test','Curve_Fitting_Toolbox');
if cfTbxLicensed
    fitCubicSpline = @(x,y) csaps(x,y);
else
    fitCubicSpline = @(x,y) unlic_csaps(x,y);
end

splinePolynomial  = fitCubicSpline(xControlPoints, yControlPoints);
padBaselineSignal = ppval(splinePolynomial, 1:padDataLenght);

baselineSignal = padBaselineSignal(extraLength:end-extraLength-1);

% de-trend signal
detrendedData = inputData - baselineSignal;