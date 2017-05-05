function meshToSTLfileASCII ( verts, faces, fileName )


fid = fopen( fileName, 'w' );
fprintf(fid,'%s\n', 'solid carlosMesh');

numOfFacets = length(faces);

for thisFacet = 1 : numOfFacets

%   Get vertex    
vertex1 = verts ( faces ( thisFacet, 1 ), : ) ;
vertex2 = verts ( faces ( thisFacet, 2 ), : ) ;
vertex3 = verts ( faces ( thisFacet, 3 ), : ) ;    

%   Compute triangle's normal
Uvector = vertex2 - vertex1 ;
Vvector = vertex3 - vertex1 ;

Nx = Uvector(2)*Vvector(3) - Uvector(3)*Vvector(2) ;
Ny = Uvector(3)*Vvector(1) - Uvector(1)*Vvector(3) ;
Nz = Uvector(1)*Vvector(2) - Uvector(2)*Vvector(1) ;

SurfaceNormal = [ Nx, Ny, Nz ];

fprintf(fid,'%s %f %f %f\n', 'facet normal', SurfaceNormal(1), SurfaceNormal(2), SurfaceNormal(3) );

fprintf(fid,'%s\n','outer loop');
fprintf(fid,'%s %f %f %f\n','vertex', vertex1(1), vertex1(2), vertex1(3));
fprintf(fid,'%s %f %f %f\n','vertex', vertex2(1), vertex2(2), vertex2(3));
fprintf(fid,'%s %f %f %f\n','vertex', vertex3(1), vertex3(2), vertex3(3));
fprintf(fid,'%s\n','endloop');

fprintf(fid,'%s\n','endfacet');

end

fprintf(fid,'%s\n','endsolid carlosMesh');
fclose(fid);