
function [videoData, videoInfo] = loadVideo ( videoFile )

%{

videoFolder = 'C:\Colours\Matlab PhD\2k15\2-Código\OpticalMappingViewer';
currentFile = '290812_NIRdye_LED_Tyrode_300fps_01.var';
videoFile   =  fullfile(videoFolder, currentFile);

%}

% helpers
validateFile = @(fileName) exist(fileName, 'file') > 0;

assert(validateFile(videoFile), 'File not found');


% save var file info
videoInfo = [];

% get code from Carolina Curiel
vFid         = fopen(videoFile, 'r');
vHeader      = fread(vFid, 6, 'uint32');
%TO-DO: check horizontal and vertical axis
vResVertical   = vHeader(2); 
vResHorizontal = vHeader(3);
numOfFrames    = vHeader(4);
framesPerSec   = vHeader(6);
vFormat        = vHeader(5)*8;
vFormatType    = sprintf('uint%d', vFormat);
vLength        = numOfFrames/framesPerSec;

% Allocate video
videoData  = zeros(vResVertical, vResHorizontal, numOfFrames);
outputSize = [vResHorizontal, vResVertical];
% TO-DO: Find a cleverer way to do this
for frameIdx = 1:numOfFrames
    videoData(:,:,frameIdx) = fread(vFid, outputSize, vFormatType)';
end

fclose(vFid);

videoInfo.vResVertical   = vResVertical;
videoInfo.vResHorizontal = vResHorizontal;
videoInfo.numOfFrames    = numOfFrames;
videoInfo.framesPerSec   = framesPerSec;
videoInfo.vLength        = vLength;

[~, fileName]      = fileparts(videoFile);
videoInfo.fileName = fileName;

% figure,
% imshow(videoData(:,:,frameIdx), [])
