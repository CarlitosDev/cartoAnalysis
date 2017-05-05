function filteredData = filterChebyShev(inputSignal, samplingFreq, varargin)

%% filterChebyShev filter signal using Chebyshev Type II
% TO-DO: Add method information
% 
% Usage
%   filteredData = filterChebyShev(inputSignal, samplingFreq, ...
%    'filterFrequency' , [0.5, 50], ...
%    'filterType'      ,'bandpass', ...
%    'filterOrder'     ,        12, ...
%    'filterRipple'    ,        20, ...
%    'showFilterDesign',      true);
% Carlos Aguilar - April 1st 2k15


% Process input arguments.
  p = inputParser;
  p.FunctionName = 'filterChebyShev';

  vDouble     = @(x) isa(x, 'double');
  vLogical    = @(x) islogical(x);
  vFilterType = @(x) find(strcmp(x, {'bandpass', 'high', 'low', 'stop'}));

  p.addRequired  ('inputSignal'     ,              vDouble);
  p.addRequired  ('samplingFreq'    ,              vDouble);
  p.addOptional  ('filterFrequency' ,  [0.2, 60],  vDouble);
  p.addOptional  ('filterType'      , 'bandpass',  vFilterType);
  p.addOptional  ('filterOrder'     ,          8,  vDouble);
  p.addOptional  ('filterRipple'    ,         20,  vDouble);
  p.addOptional  ('showFilterDesign',      false,  vLogical);

  p.parse(inputSignal, samplingFreq, varargin{:});

  inputSignal      = p.Results.inputSignal;
  samplingFreq     = p.Results.samplingFreq;
  filterFrequency  = p.Results.filterFrequency;
  filterType       = p.Results.filterType;
  filterOrder      = p.Results.filterOrder;
  filterRipple     = p.Results.filterRipple;
  showFilterDesign = p.Results.showFilterDesign;


%% set up filter

% Adjust order if bandpass or stopband (cheby2 returns an order 2*n
% bandpass analog filter )
  vFilterDesign = ...
    @(filterFrequency, filterType) numel(filterFrequency)==2 & any(strcmp(filterType, {'bandpass', 'stop'}));

  if vFilterDesign(filterFrequency, filterType)
    actualOrder = filterOrder/2;
    freqText    = sprintf('(%3.2f,%3.2f)', filterFrequency(1), filterFrequency(2));
  else
    freqText    = sprintf('%3.2f', filterFrequency);
    actualOrder = filterOrder;
  end

  Wn = filterFrequency/(samplingFreq/2);
  
  [numSignals, numSamples] = size(inputSignal);
  
  assert(numSignals<numSamples, 'Arrange data in (signal x samples) format\n');

%% Cheby filter 

  % malloc
  filteredData = zeros(numSignals, numSamples);


  % (i) Get poles, zeros and gain
  %     % get filter $H(z) = fract{B(z)}{A(z)}$ being B(z) poles
  %     [fZeros, fPoles, fGain] = cheby2(actualOrder, filterRipple, Wn, filterType);
  %     [sos,g] = zp2sos(fZeros, fPoles, fGain);% Convert to SOS form
  %     [b,a]   = zp2tf(fZeros, fPoles, fGain);
  
  % (ii) Easier: get transfer function
  [filtNum, filtDenom] = cheby2(actualOrder, filterRipple, Wn, filterType);
  
  for idx = 1:numSignals
    filteredData(idx, :) = filtfilt(filtNum, filtDenom, inputSignal(idx, :));
  end
  
  % Show filter design (from Matlab help)
  if showFilterDesign
    [sos,g] = tf2sos(filtNum, filtDenom);
    Hd      = dfilt.df2tsos(sos,g); % Create a dfilt object
    h       = fvtool(Hd);           % Plot magnitude response
    set(h,'Analysis','freq')        % Display frequency respon
  end

  fprintf('\tfilterChebyShev: Filtering %dx%d signals with \n\tChevTII order %d, att %d, type %s %s Hz\n', ...
    numSignals, numSamples, filterOrder, filterRipple, filterType, freqText);
  
% TO-DO: Add here code for svn/gitHub  
% $Revision: $
% $Author: carlosAguilar $
% $Date: $
% Copyright 2015 OGTel LTD.