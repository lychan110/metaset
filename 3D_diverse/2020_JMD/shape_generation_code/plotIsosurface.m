function fv = plotIsosurface(geometry)
% ==============================================================================
% Plot logical voxel or signed distance function data as isosurface
% 
% Input: geometry - three options:
%                   1) 3D logical matrix of voxels;
%                   2) 3D matrix of signed distance field values (aka level-set), 
%                      where phi > 0 is solid material;
%                   3) struct of geometry.vertices and 'geometry.faces (both nx3 array).
%
% Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020.
% "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“
% Journal of Mechanical Design. March 2021; 143(3): 031707.
% 
% Author: Yu-Chin Chan (ychan@u.northwestern.edu), 5/2/2019
% Late updated: 5/11/2019, 6/18/2019
% ==============================================================================
if isstruct(geometry)
    fv = geometry;
else
    if islogical(geometry)
        fv = isosurface(padarray(geometry, [1 1 1]), 0.1);
    else
        fv = phi2fv(geometry, 1, 1, 0);
    end
end
p = patch(fv);
p.FaceColor = 'red';
p.EdgeColor = 'none';
xlabel('x'), ylabel('y'), zlabel('z')
set(gca, 'fontsize', 14)
view(3);
axis equal
axis tight
camlight 
lighting gouraud
end