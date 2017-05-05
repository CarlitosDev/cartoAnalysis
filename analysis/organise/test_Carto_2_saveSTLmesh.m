% %% Run fastICA over a collection of signals
% 
% This script runs an independent component analysis -fastICA- on ECG/EGM
% data given a LISFTOFSIGNALS and a LISTOFPOINTS.
% 
% To set up the script, follow the section below.
% 
%
%%

listOfSelectedIcs = {[1,2]};
%listOfSelectedIcs = {1, 2, 3, [1,2], [1,2,3]};

testType = 'NewDataS2';

% carto file (also relative to main folder)
%     path2CartoFile     = fullfile('data','cartoXP','handles.mat');   
%     fullPath2CartoFile = fullfile(mainFolder, path2CartoFile);
%path2CartoFile = fullfile('data','Auriculas','#4', 'Carto2', 'EG.mat');
%/Users/carlosAguilar/Documents/Matlab PhD/codigoCarlos/data/Auriculas/#5/Carto/1-Map
%path2CartoFile = fullfile('data','Auriculas','#5', 'Carto', '1-Map', 'EG.mat');
%path2CartoFile = fullfile('data','Auriculas','#4', 'Carto1', 'EG.mat');
path2CartoFile = fullfile('data','Auriculas','#2', 'Carto', 'EG.mat');
fullPath2CartoFile = fullfile(rootFolder(), path2CartoFile); 
        

%% load CartoXP

[data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile);


%%