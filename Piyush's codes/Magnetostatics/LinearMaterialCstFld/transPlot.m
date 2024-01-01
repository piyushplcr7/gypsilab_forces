% Custom plotting function for interior and exterior meshes
function [] = transPlot(bndmesh)
    axis equal;
%     H = patch('Faces',bndmesh_e.elt, 'Vertices',bndmesh_e.vtx,'Marker','o');
    H = patch('Faces',bndmesh.elt, 'Vertices',bndmesh.vtx);
    colorData = bndmesh.col;
    set(H,'FaceVertexCData',colorData,'FaceColor','flat','FaceAlpha',0.2);
%     hold on;
%     plot(bndmesh_i);
end