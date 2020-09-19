function [mesh1,mesh2] = mshSplit(mesh,X0,U,varargin)
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
%|    #    |   FILE       : mshSplit.m                                    |
%|    #    |   VERSION    : 0.51                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2019                                    |
%|  / 0 \  |   LAST MODIF : 21.06.2019                                    |
%| ( === ) |   SYNOPSIS   : Split mesh with planar cut                    |
%|  `---'  |                (to be clearly improved)                      |
%+========================================================================+

p = inputParser;
p.addOptional('stable',true)
% Stable true : don't modify elements, split according to element center
% positions.
% Stable false : make a clean cut. Will lead to adding new vertices to
% resolve the cut. 

p.parse(varargin{:});
stab = p.Results.stable;

N = U./norm(U);
height = @(list_X)((list_X - ones(size(list_X,1),1)*X0)*N');
is_above = @(list_X)( height(list_X)> 0);
        % true : above the plane. false : below the plane
is_below = @(list_X)(~is_above(list_X));
        
if stab || strcmp(mesh.type,'point')
    % Normalize direction
    
    % Split mesh
    I = is_above(mesh.ctr);
    mesh1  = mesh.sub(I);
    mesh2  = mesh.sub(~I);
else
    switch mesh.type
        case 'point'
            % Done above.
        case 'segment'
            [A,B] = ABCD(mesh);
            % Find elements that are entirely above
            above = and(is_above(A),is_above(B));
            % Find elements that are entirely below
            below = and(is_below(A),is_below(B));
            % Find element that are in between
            in_between = ~ or(above,below);
            % Number of elements in between
            nint = sum(in_between);
            
            Aint = A(in_between,:); Bint = B(in_between,:);
            vtx_int = intersectSegPlane(Aint,Bint,X0,N); % Find intersection
            vtx = [mesh.vtx;vtx_int]; % we add the intersection points to the vertices
            elt = mesh.elt(or(above,below),:); % remove elements in between
            col = mesh.col(or(above,below),:);
            elt = [elt; mesh.elt(in_between,1), mesh.nvtx + (1:nint)']; % add their split versions
            elt = [elt; mesh.nvtx + (1:nint)', mesh.elt(in_between,2)];
            col = [col; mesh.col(in_between); mesh.col(in_between)]; % inherit colors
            tmp = msh(vtx,elt,col); 
            % Once the cutting plane is resolved, we can use the stable
            % option, because no element is "across" the plane. 
            [mesh1,mesh2] = mshSplit(tmp,X0,U,'stable',true);
        case 'triangle'
            [A,B,C] = ABCD(mesh);
            above = all([is_above(A),is_above(B),is_above(C)],2);
            below = all([is_below(A),is_below(B),is_below(C)],2);
            in_between = find(~or(above,below));
            nint = length(in_between);
            indIJK = mesh.elt(in_between,:);
            I = zeros(nint,3); J = I; K = I;
            indi = zeros(nint,1); indj = indi; indk = indi;
            indI = indi; indJ = indi; indK = indi;
            for i = 1:nint
                ABC_i = [A(in_between(i),:);B(in_between(i),:);C(in_between(i),:)];
                above_i = is_above(ABC_i);
                if sum(above_i)==1
                    [~,indi(i)] = max(above_i);
                    I(i,:) = ABC_i(indi(i),:);
                    indj(i) = mod(indi(i),3)+1; indk(i) = mod(indj(i),3) + 1;
                    J(i,:) = ABC_i(indj(i),:);
                    K(i,:) = ABC_i(indk(i),:);
                elseif sum(above_i) ==2
                    [~,indi(i)] = min(above_i);
                    I(i,:) = ABC_i(indi(i),:);
                    indj(i) = mod(indi(i),3)+1; indk(i) = mod(indj(i),3) + 1;
                    J(i,:) = ABC_i(indj(i),:);
                    K(i,:) = ABC_i(indk(i),:);
                end
                indI(i) = indIJK(i,indi(i));
                indJ(i) = indIJK(i,indj(i));
                indK(i) = indIJK(i,indk(i));
                % Intersection points are located on IJ and IK.
            end
            vtx_IJ = intersectSegPlane(I,J,X0,N);
            vtx_IK = intersectSegPlane(I,K,X0,N);
            vtx = [mesh.vtx;vtx_IJ;vtx_IK];
            indIJ = mesh.nvtx + (1:nint)'; 
            indIK = mesh.nvtx + nint + (1:nint)';
            elt = mesh.elt(or(above,below),:);
            elt = [elt; [indI, indIJ, indIK]]; 
            elt = [elt; [indJ,indIK,indIJ]];
            elt = [elt; [indJ,indK,indIK]];
            col = mesh.col(or(above,below),:);
            col = [col; repmat(mesh.col(in_between),4,1)];
            tmp = msh(vtx,elt,col);
            [mesh1,mesh2] = mshSplit(tmp,X0,N,'stable',true);
        case 'tetrahedron'
            error('Not available yet. Feel free to implement it !')
    end
end



end

function[x] = intersectSegPlane(A,B,X0,N)
    u = sum((A - ones(size(A,1),1)*X0)*((ones(size(A,1),1)*N)'),2)./...
                sum((A - B)*((ones(size(A,1),1)*N)'),2);
    x = A + u*ones(1,3).*(B - A);
end

