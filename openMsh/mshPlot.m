function mshPlot(mesh,color)
%+========================================================================+
%|                                                                        |
%|                 OPENMSH - LIBRARY FOR MESH MANAGEMENT                  |
%|           openMsh is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2017-2018.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : mshPlot.m                                     |
%|    #    |   VERSION    : 0.40                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 14.03.2018                                    |
%| ( === ) |   SYNOPSIS   : Plot mesh and data                            |
%|  `---'  |                                                              |
%+========================================================================+


if ~exist('color','var')
    color = 'b';
end
if isnumeric(color)&&isscalar(color)
    color = repmat(color,mesh.nvtx,1);
end

% Patch
H = patch('Faces',mesh.elt, 'Vertices',mesh.vtx,'Marker','o');

switch mesh.type
    case 'point'
        if ischar(color)
            set(H,'MarkerFaceColor',color);
        else
            set(H,'MarkerFaceColor','flat','FaceVertexCData',color);
        end
    case 'segment'
        if ischar(color)
            set(H,'MarkerFaceColor',color);
            set(H,'MarkerEdgeColor','k');
            set(H,'EdgeColor',color);
        else
            set(H,'MarkerFaceColor','flat','FaceVertexCData',color);
            set(H,'EdgeColor','interp');
            set(H,'MarkerEdgeColor','k');
        end
        
    case 'triangle'
        if ischar(color)
            set(H,'Marker','.','MarkerFaceColor','k','EdgeColor','k','FaceColor',color)
        else
            set(H,'Marker','none','EdgeColor','none','FaceColor','interp','FaceVertexCData',color)
        end
        
    case 'tetrahedron'
        
end

return

end
