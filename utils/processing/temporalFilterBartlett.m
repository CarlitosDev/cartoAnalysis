
function videoDataFiltered = temporalFilterBartlett( videoData, videoInfo, kernelSize ) 
%{

videoFolder = 'C:\Colours\Matlab PhD\2k15\2-Código\OpticalMappingViewer';
currentFile = '290812_NIRdye_LED_Tyrode_300fps_01.var';
videoFile   =  fullfile(videoFolder, currentFile);

%}

%% Get Bartlett mask

% $ w(i,j) = 1 - \sqrt{\fract{k^2+l^2}{m^2}} $

% Code from Carolina Curiel
% Máscara cónica de Bartlett (Mironov 2006)
n_mask = kernelSize;
mask   = zeros(n_mask, n_mask);

for i=1:n_mask
    for j=1:n_mask
        k=floor(i-n_mask/2);
        l=floor(j-n_mask/2);
        m=n_mask;
        resta=min(sqrt((k^2+l^2)/m^2),1);
        mask(i,j)=1-resta;
    end
end

%% filter video
%videoDataFiltered = imfilter( videoData, mask, 'symmetric');
% Read the video as a signal over time

videoDataFiltered = zeros(size(videoData));


NFFT = 2^nextpow2(videoInfo.numOfFrames);

maskFrequency = fft(mask(2,:), NFFT)';

for i=1:videoInfo.vResHorizontal
  for j=1:videoInfo.vResVertical
        
        profileFFT = fft(squeeze(videoData(i,j,:)), NFFT);
        filteredProfile = ifft(profileFFT.*maskFrequency);
        videoDataFiltered(i, j, :) = filteredProfile(1:videoInfo.numOfFrames) - mean(filteredProfile(1:videoInfo.numOfFrames));
    
   end
end

% figure,
% plot(squeeze(videoData(i,j,:)))
% hold on,
% plot(squeeze(videoDataFiltered(i,j,:)), 'r')

% %% 
% for i=1:vResHorizontal
%   for j=1:vResVertical
%         
%         filteredProfile = conv(squeeze(videoData(i,j,:)), mask(2,:), 'same');
%     
%    end
% end