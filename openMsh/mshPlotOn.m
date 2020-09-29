function[] = mshPlotOn(mesh,func)
% Plot a function defined on the mesh.
% Patch



if isa(func,'function_handle')
    colorData = func(mesh.vtx);
else
    assert(size(func,1) == mesh.nvtx);
    colorData = func;
end


if is1d(mesh)
    plot(mesh); hold on
    mesh.vtx(:,2) = colorData;
    plot(mesh.vtx(:,1),mesh.vtx(:,2));
    return
% elseif is2d(mesh)
%     if strcmp(mesh.type,'segment')
%         plot(mesh); hold on
%     end
%     view(4,14);
%     mesh.vtx(:,3) = colorData;  
% end
end

H = patch('Faces',mesh.elt, 'Vertices',mesh.vtx,'Marker','o');



switch mesh.type
    case 'point'
        if ischar(colorData)
            set(H,'MarkerFaceColor',colorData);
        else
            set(H,'MarkerFaceColor','flat','FaceVertexCData',colorData);
        end
    case 'segment'
        if ischar(colorData)
            set(H,'Marker','x','MarkerEdgeColor','k','MarkerSize',5);
            set(H,'EdgeColor','interp','FaceVertexCData',colorData);
        else
            set(H,'EdgeColor','flat','FaceVertexCData',colorData,'LineWidth',2);
            set(H,'Marker','none');
            
        end
        
    case 'triangle'
        if ischar(colorData)
            set(H,'Marker','.','MarkerFaceColor','k','FaceColor',colorData)
        else
            set(H,'FaceVertexCData',colorData,'FaceColor','interp');
            set(H,'EdgeColor','none');
            set(H,'Marker','.');
        end
        
    case 'tetrahedron'
        
        m = mesh.bnd; % as we cannot se the interior...
        
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