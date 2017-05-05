%% my code

addpath('D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\carlosUtils');
addpath('D:\Work\MATLABPromotionsModel\Matlab PhD\matlab code\toolboxes\FastICA_2.5');
addpath(genpath('D:\Work\MATLABPromotionsModel\mDev\mMatlab'));

mMeshFolder = 'D:\Work\MATLABPromotionsModel\Matlab PhD\TesisCarlos\3-CodigoIcaFA - Aitor\3-CodigoIcaFA';

% that loads 'handles'
load('D:\Work\MATLABPromotionsModel\Matlab PhD\TesisCarlos\3-CodigoIcaFA - Aitor\3-CodigoIcaFA\handles.mat')
fieldnames(handles)


information = handles.DB;
data        = information.EG;


meshVertex = [meshData.X; meshData.Y; meshData.Z]';
meshFaces  = meshData.faces;

figure,
h = trisurf(meshFaces, ...
    meshVertex(:, 1), meshVertex(:, 2), meshVertex(:, 3));
pointsInfoStr = num2str(pointsInfo');
text(meshVertex(:, 1), meshVertex(:, 2), meshVertex(:, 3), ...
    pointsInfoStr);

get(h)



% save new one
fileName = fullfile(mMeshFolder,'leftAtria.stl');
tic
SaveSTLfileASCII ( meshFaces, meshVertex, fileName );
toc

tic
[mfRefinedMesh, mnTriangulation] = LoopSubdivision( meshVertex, meshFaces );
fprintf (mFid, '\nDone in %3.f seconds',toc);
% save new one
fileName = [mMeshFolder,'\', mMesh1(1:end-4), 'Refined.stl' ];
tic
SaveSTLfileASCII ( mnTriangulation, mfRefinedMesh, fileName );
toc



figure

trisurf(meshFaces, meshVertex(:, 1), meshVertex(:, 2), ...
    meshVertex(:, 3), meshVertex(:, 4))