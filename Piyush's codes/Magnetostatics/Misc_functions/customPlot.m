% Custom plotting function for interior and exterior meshes
function [] = customPlot(bndmesh_i,bndmesh_e)
    axis equal;
%     H = patch('Faces',bndmesh_e.elt, 'Vertices',bndmesh_e.vtx,'Marker','o');
    H = patch('Faces',bndmesh_e.elt, 'Vertices',bndmesh_e.vtx);
    colorData = bndmesh_e.col;
    set(H,'FaceVertexCData',colorData,'FaceColor','flat','FaceAlpha',0.2);
    hold on;
    plot(bndmesh_i);
end