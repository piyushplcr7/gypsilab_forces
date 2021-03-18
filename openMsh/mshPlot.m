function mshPlot(mesh,colorData)

if ~exist('colorData','var')||isempty(colorData)
    colorData = mesh.col;
end

if isa(colorData,'function_handle')
    colorData = colorData(mesh.vtx);
end


assert(or(size(colorData,1) == length(mesh),ischar(colorData)));
H = patch('Faces',mesh.elt, 'Vertices',mesh.vtx,'Marker','o');

axis equal
switch mesh.type
    case 'point'
        if ischar(colorData)
            set(H,'MarkerFaceColor',colorData);
        else
            set(H,'MarkerFaceColor','flat','FaceVertexCData',colorData);
        end
    case 'segment'
        if ischar(colorData)
            set(H,'Marker','x','MarkerEdgeColor',colorData,'MarkerSize',5);
            set(H,'EdgeColor',colorData,'LineWidth',2);
        else
            color2 = colorData;
            set(H,'EdgeColor','flat','FaceVertexCData',color2,'LineWidth',2);
            set(H,'Marker','x','MarkerEdgeColor','k','MarkerFaceColor','none','MarkerSize',5);
        end
        set(H,'Marker','none');
        
    case 'triangle'
        if ischar(colorData)
            set(H,'Marker','.','MarkerFaceColor','k','EdgeColor','k','FaceColor',colorData)
        else
            set(H,'FaceVertexCData',colorData,'FaceColor','flat');
            set(H,'Marker','.');
        end
        
    case 'tetrahedron'
        
        m = mesh.bnd; % as we cannot see the interior...
        
        for i = 1:mesh.nelt % In this loop : assign to each face the color
            % of one associated tetrahedron that it touches. Faces will be
            % coloured several times but at least every one gets a color in
            % the end.
            if ischar(colorData)
                delete(H);
                plot(m,colorData);
            else
                plot(m);
            end
        end
end
end