function [detrendedData, baselineSignal]=removeBaselineWandering(inputData)

% REMOVEBASELINEWANDERING Remove wandering of the base line of an input
% signal.
%
% 

numberOfSamples = numel(inputData);
% numSmpFactor    = factor(numberOfSamples);
% [uniqueFactors, ~, ic] = unique(numSmpFactor);
% add a snippet to do it automatically.
% So far, let's pick 5^3

stepSize = 5^3;
xControlPoints    = 0:stepSize:numberOfSamples;
xControlPoints(1) = 1;

% let's smooth the signal before getting the values for the control points
filterKernel = ones(1, stepSize)./stepSize;
filteredInputData = conv(inputData, filterKernel, 'same');
yControlPoints    = filteredInputData(xControlPoints);



% fit the control point to a cubic spline
% work this out...
cfTbxLicensed = license('test','Curve_Fitting_Toolbox');
if cfTbxLicensed
    fitCubicSpline = @(x,y) csaps(x,y);
else
    fitCubicSpline = @(x,y) unlic_csaps(x,y);
end

splinePolynomial = fitCubicSpline(xControlPoints, yControlPoints);
baselineSignal   = ppval(splinePolynomial, 1:numberOfSamples);

% de-trend signal
detrendedData = inputData - baselineSignal;