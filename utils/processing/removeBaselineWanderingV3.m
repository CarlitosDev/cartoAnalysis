function [detrendedData, baselineSignal] = removeBaselineWanderingV3(inputData, samplingPercentage)

% REMOVEBASELINEWANDERING Remove wandering of the base line of an input
% signal using a spline approximation.
%
% This version -unlike the previous ones - doesn't depend on the curve fitting toolbox 
%
%   Carlos Aguilar Nov 2k16

%%
  numberOfSamples = numel(inputData);
  extraLength     = round(numberOfSamples*0.2);
  extraLength     = 0.5*(extraLength + double(mod(extraLength,2)~=0));

  % symmetric padding (padarray belongs to image processing tbx)
  leftSide      = fliplr(inputData(1:extraLength));
  rightSide     = fliplr(inputData(end-extraLength+1:end));
  padInputData  = horzcat(leftSide, inputData, rightSide);
  
  padDataLenght = numel(padInputData);

  % Set the step to some percent of the data
  stepSize = floor(samplingPercentage*padDataLenght/100);

  % add a snippet to do it automatically.
  xControlPoints    = 0:stepSize:padDataLenght;
  xControlPoints(1) = 1;

  stepsVector  = 1:padDataLenght;


  % let's smooth the signal before getting the values for the control points
  filterKernel      = ones(1, stepSize)./stepSize;
  filteredInputData = conv(padInputData, filterKernel, 'same');
  yControlPoints    = filteredInputData(xControlPoints);

  splineInterp      = spline(stepsVector    , filteredInputData, xControlPoints);
  padBaselineSignal = interp1(xControlPoints, splineInterp, stepsVector, 'spline');

  baselineSignal = padBaselineSignal(extraLength:end-extraLength-1);

  % de-trend signal
  detrendedData = inputData - baselineSignal;