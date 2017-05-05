%% readDataAndSave

mainFolder     =  fullfile(rootFolder, 'data','cartoXPcarlos');
path2CartoFile = fullfile('#2', 'Carto', 'EG.mat');
path2CartoFile = fullfile('#3', 'Carto', 'EG.mat');
path2CartoFile = fullfile('#4', 'Carto1', 'EG.mat');
path2CartoFile = fullfile('#4', 'Carto2', 'EG.mat');
path2CartoFile = fullfile('#5', 'Carto', 'EG.mat');


fullPath2CartoFile = fullfile(mainFolder, path2CartoFile);


[pathStr, ~, ~] = fileparts(fullPath2CartoFile);

%%
[data, pointsInfo, meshData, cartoData] = loadCartoXPData(fullPath2CartoFile);
%save(fullfile(pathStr, 'cartoData.mat'), 'data', 'pointsInfo', 'meshData', 'cartoData');

%% save STL file
% stlFileName = fullfile(pathStr, 'cartoDataMesh.stl');
% meshToSTLfileASCII( meshData.meshVertex, meshData.faces, stlFileName );

%% plot mesh and carto points

% arrange data
xData = cartoData.pointsPosition(:,1);
yData = cartoData.pointsPosition(:,2);
zData = cartoData.pointsPosition(:,3);

figPathName  = fullfile(pathStr, 'meshPlot.jpg');
currentFaces = meshData.faces;
currentVert  = [meshData.X,meshData.Y,meshData.Z];

%%
figure,
% plot carto points
plot3(xData, yData, zData, ...
    'o', 'MarkerSize' , 12, ...
    'MarkerFaceColor',[.49 1 .63]);

% plot mesh
hold on,
trimesh(currentFaces, ...
        meshData.X, meshData.Y, meshData.Z, ...
        'FaceColor', 'none', ...
        'EdgeColor', [0 0 1]);

figPathName  = fullfile(pathStr, 'meshPlot.jpg');

set(gcf, 'Position', [2649, 189, 704, 808]);
    
hax = gca();
figTitle = strrep(path2CartoFile, '\', ' ');
title(hax, figTitle);  
% hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'jpeg');
% close(gcf);


%%
% let's see if the point-facet thing makes any sense
figPathName  = fullfile(pathStr, 'meshPlotPoints.jpg');

figure,
% plot carto points
plot3(xData, yData, zData, ...
    'o', 'MarkerSize' , 12, ...
    'MarkerFaceColor',[.49 1 .63]);


% plot mesh
hold on,
patch('Faces'          , currentFaces         , ...
      'Vertices'       , currentVert          , ...
      'FaceVertexCData', meshData.facetColours, ...
      'FaceColor'      , 'flat'               , ...
      'facealpha'      , 0.8                  , ...      
      'CDataMapping'   , 'direct');
          
set(gcf, 'Position', [2649, 189, 704, 808]);

hax = gca();
figTitle = strrep(path2CartoFile, '\', ' ');
title(hax, figTitle);  
% hgexport(gcf, figPathName, hgexport('factorystyle'), 'Format', 'jpeg');
% close(gcf);

    
 