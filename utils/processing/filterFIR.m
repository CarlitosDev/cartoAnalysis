function filteredData = filterFIR(inputSignal, samplingFreq, varargin)

%% filterFIR filter signal using FIR 
% TO-DO: Add method information
% 
% Usage
%   filteredData = filterFIR(inputSignal, samplingFreq, ...
%    'filterFrequency' , [0.5, 50], ...
%    'filterType'      ,'bandpass', ...
%    'filterOrder'     ,        12, ...
%    'filterRipple'    ,        20, ...
%    'showFilterDesign',      true);
% Carlos Aguilar - June 12th 2k15


% Process input arguments.
  p = inputParser;
  p.FunctionName = 'filterFIR';

  vDouble     = @(x) isa(x, 'double');
  vLogical    = @(x) islogical(x);
  vFilterType = @(x) find(strcmp(x, {'bandpass', 'high', 'low', 'stop'}));

  p.addRequired  ('inputSignal'     ,              vDouble);
  p.addRequired  ('samplingFreq'    ,              vDouble);
  p.addOptional  ('filterFrequency' ,  [0.2, 60],  vDouble);
  p.addOptional  ('filterType'      , 'bandpass',  vFilterType);
  p.addOptional  ('filterOrder'     ,        129,  vDouble);
  p.addOptional  ('filterRipple'    ,         20,  vDouble);
  p.addOptional  ('showFilterDesign',      false,  vLogical);

  p.parse(inputSignal, samplingFreq, varargin{:});

  inputSignal      = p.Results.inputSignal;
  samplingFreq     = p.Results.samplingFreq;
  filterFrequency  = p.Results.filterFrequency;
  filterType       = p.Results.filterType;
  filterOrder      = p.Results.filterOrder;
  showFilterDesign = p.Results.showFilterDesign;


%% set up filter

% Adjust order if bandpass or stopband
  vFilterDesign = ...
    @(filterFrequency, filterType) numel(filterFrequency)==2 & any(strcmp(filterType, {'bandpass', 'stop'}));

  %makeOrderOdd = @(x) x+double(mod(x,2) == 0);
  makeOrderEven = @(x) x+double(mod(x,2) == 1);

  if vFilterDesign(filterFrequency, filterType)
    freqText = sprintf('(%3.2f,%3.2f)', filterFrequency(1), filterFrequency(2));
  else
    freqText = sprintf('%3.2f', filterFrequency);
  end
  
  actualOrder = makeOrderEven(filterOrder);
  
  Wn = filterFrequency/(samplingFreq/2);
  
  [numSignals, numSamples] = size(inputSignal);
  
  assert(numSignals<numSamples, 'Arrange data in (signal x samples) format\n');

%% FIR filter 

  % malloc
  filteredData = zeros(numSignals, numSamples);

  
  % Get transfer function
  filtDenom = 1;
  filtNum   = fir1(actualOrder, Wn, filterType);

  
  for idx = 1:numSignals
    filteredData(idx, :) = filtfilt(filtNum, filtDenom, inputSignal(idx, :));
  end
  
  % Show filter design (from Matlab help)
  if showFilterDesign
    [sos,g] = tf2sos(filtNum, filtDenom);
    Hd      = dfilt.df2tsos(sos,g); % Create a dfilt object
    h       = fvtool(Hd);           % Plot magnitude response
    set(h,'Analysis','freq')        % Display frequency response
  end

  fprintf('\tfilterFIR: Filtering %dx%d signals with \n\tFIR order %d, type %s %s Hz\n', ...
    numSignals, numSamples, filterOrder, filterType, freqText);
  
% TO-DO: Add here code for svn/gitHub  
% $Revision: $
% $Author: callosAguilar $
% $Date: $
% Copyright 2015 OGTel LTD.