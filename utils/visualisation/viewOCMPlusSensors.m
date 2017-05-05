function viewOCMPlusSensors ( videoData, videoInfo, sensorsXY )

% TO-DO: Add some basic checks here for matrix dimensions
% MAke sure sensorsXY is in the right format

updateImage = 1/3;% secs

%{

sensorsXY = [32,32 ; 30, 30; 28, 28];

%}


% Set sensors
[numSensors, ~] = size(sensorsXY);
sensorIntensity = zeros(videoInfo.numOfFrames, numSensors);

for idx=1:numSensors
    sensorIntensity(:, idx) = squeeze( videoData(sensorsXY(idx,1), sensorsXY(idx,2),:));
end




%% Plot the thing
figure,
h1 = subplot(121);
% frameIdx = 1;
% imshow(videoData(:,:,frameIdx), [])
numSensor = 1;
h2 = subplot(122);
plot(h2, sensorIntensity(:, numSensor), 'Color', [1 1 1]);
hold on,
for frameIdx = 1:videoInfo.numOfFrames
    subplot(121), imshow(videoData(:,:,frameIdx), [])
    plot(h2, 1:frameIdx, sensorIntensity(1:frameIdx, numSensor), 'Color', [0 0 1]);
    hold on
    plot(h2, 1:frameIdx, sensorIntensity(1:frameIdx, numSensor+1), 'Color', [0 1 0]);
    hold on
   pause(updateImage)
end