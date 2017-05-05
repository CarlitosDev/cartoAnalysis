%% initialise
addpath(genpath(rootFolder()));

%% add the folders left out from GIT
% external code from Marga & JL
extCode = fullfile(rootFolder(),'..','code (ungitted)');
if isdir(extCode)
    addpath(genpath(extCode));
end

% private data
dataFolder = fullfile(rootFolder(),'..','Data (ungitted)');
if isdir(dataFolder)
    addpath(genpath(dataFolder));
end

%% edit the most recent ...

cd(rootFolder());