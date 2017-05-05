function addWorkDirectory2Path()


    mainFolder  = rootFolder();
    % add utils to path
    currentFolder = 'utils';
    addpath(fullfile(mainFolder,currentFolder));

    % add fastICA     
    fastICAFolder = fullfile('external code', 'FastICA_2.5');
    addNewPath(fullfile(mainFolder, fastICAFolder));


    % move to tests folder
    currentFolder = 'tests';
    cd(fullfile(mainFolder, currentFolder));