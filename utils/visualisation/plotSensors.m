function [sensorProfile, figAxes] = plotSensors ( videoData, videoInfo, sensorsXY )

% TO-DO: Add some basic checks here for matrix dimensions
% Make sure sensorsXY is in the right format and is less than 64 elements
% (use mod otherwise to repeat colour)


%{
    sensorsXY = [32,32 ; 30, 30; 28, 28];
    [sensorProfile, figAxes] = plotSensors ( videoData, videoInfo, sensorsXY )
%}

% Set graphics here
plotColours = colormap('Lines');
h = figure();
figAxes = axes();

% Set sensors
[numSensors, ~] = size(sensorsXY);
sensorProfile = zeros(videoInfo.numOfFrames, numSensors);


% Get and plot intensity profiles
for idx = 1:numSensors
    
    sensorProfile(:, idx) = squeeze( videoData(sensorsXY(idx,1), sensorsXY(idx,2),:));
    % plot    
    plot(figAxes, sensorProfile(:, idx), 'Color', plotColours(idx, :) );
    hold on;
end