
function videoDataFiltered = spatialFilterBartlett( videoData, videoInfo, kernelSize ) 
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


% TO-DO
% review this...Bartlett should be linear. A 3-by-3 should look like
% a = [0.5,1,0.5]
% b = a'*a;
% mask = b./max(b(:));

%% filter video

videoDataFiltered = imfilter( videoData, mask, 'symmetric');