%% my code

mMeshFolder = 'D:\Work\MATLABPromotionsModel\Matlab PhD\TesisCarlos\3-CodigoIcaFA - Aitor\3-CodigoIcaFA';

% that loads 'handles'
load('D:\Work\MATLABPromotionsModel\Matlab PhD\TesisCarlos\3-CodigoIcaFA - Aitor\3-CodigoIcaFA\handles.mat')
fieldnames(handles)


information = handles.DB;
data        = information.EG;


meshVertex = information.mesh_vertex;
meshFaces  = information.mesh_faces;


figure
trisurf(meshFaces, meshVertex(:,1), meshVertex(:,2), meshVertex(:,3), meshVertex(:,4)); % View torso

trisurf(Tri,X,Y,Z)