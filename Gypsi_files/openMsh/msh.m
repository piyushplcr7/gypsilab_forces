classdef msh
    
    properties
        vtx = [];      % VERTEX COORDINATES (3 dimensions)
        elt = [];      % ELEMENTS LIST (particles, edges, triangles or tetrahedron)
        col = [];      % ELEMENT COLOR (group number)
    end
    
    methods
        %% Class constructor %%
        
        function this = msh(varargin)
            % Read file
            if (length(varargin) == 1) && ischar(varargin{1})
                this = msh;
                file = varargin{1};
                ext  = file(end-2:end);
                if strcmp(ext,'msh')
                    [this.vtx,this.elt,this.col] = mshReadMsh(file);
                elseif strcmp(ext,'ply')
                    [this.vtx,this.elt] = mshReadPly(file);
                elseif strcmp(ext,'stl')
                    [this.vtx,this.elt] = mshReadStl(file);
                elseif strcmp(ext,'vtk')
                    [this.vtx,this.elt] = mshReadVtk(file);
                elseif strcmp(file(end-3:end),'mesh')
                    [this.vtx,this.elt,this.col] = mshReadMesh(file);
                else
                    error('msh.m : unavailable case')
                end
                if isempty(this.col)
                    this.col = zeros(size(this.elt,1),1);
                end
                
                % Inout only vertex
            elseif (length(varargin) == 1) && isnumeric(varargin{1})
                this.vtx = varargin{1};
                this.elt = (1:size(this.vtx,1))';
                this.col = zeros(size(this.elt,1),1);
                
                % Input vertex and elements
            elseif (length(varargin) == 2)
                this.vtx = varargin{1};
                this.elt = varargin{2};
                this.col = zeros(size(this.elt,1),1);
                
                % Input vertex, elements and colours
            elseif (length(varargin) == 3)
                this.vtx = varargin{1};
                this.elt = varargin{2};
                this.col = varargin{3}(:); % make sure it's a column.
            end
            
            % Clean mesh
%             this = clean(this);
        end
        
        % Clean: Remove duplicate vertices, remove unused vertices and relabel
        % elements accordingly.
        function this = clean(this,dst)
            if ~exist('dst','var')||isempty(dst)
                dst = [];
            end
            this = mshClean(this,dst);
        end
        
        
        %% Display %%
        
        function disp(this)
            space = '   ';
            fprintf('%s %s mesh with %d elements and %d vertices: \n\n',' ',this.type,this.nelt,this.nvtx);
            fprintf('%s vtx: [%dx%d double] \n',space,size(this.vtx,1),size(this.vtx,2))
            fprintf('%s elt: [%dx%d double] \n',space,size(this.elt,1),size(this.elt,2))
            fprintf('%s col: [%dx%d double] \n\r',space,size(this.col,1),size(this.col,2))
        end
        
        function plot(varargin)
            mshPlot(varargin{:});
        end
        
        function plotNrm(varargin)
            
            mesh = varargin{1};
            assert(mesh.dim <= 3);
            spc  = 'r';
            
            if (nargin == 2)
                spc = varargin{2};
            end
            Xctr = mesh.ctr;
            Vnrm = mesh.nrm;
            Vnrm = Vnrm.*(mesh.ndv*([1 1 1]));
            quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),...
                'Color',spc);
        end
        
        function plotTgt(varargin)
            mesh = varargin{1};
            spc  = 'r';
            if (nargin == 2)
                spc = varargin{2};
            end
            switch mesh.type
                case 'segment'
                    [A,B] = ABCD(mesh);
                    Xstart = A + (B-A)/3;
                    Vtgt = mesh.tgt;
                    Vtgt = Vtgt.*(mesh.ndv*([1 1 1]));
                    quiver3(Xstart(:,1),Xstart(:,2),Xstart(:,3),Vtgt(:,1),Vtgt(:,2),Vtgt(:,3),...
                        'Color',spc,'LineWidth',2);
                otherwise
                    % Orientation will be arbitrary, since an edge can be
                    % shared by two triangles...
                    edg = mesh.edg;
                    plotTgt(edg,spc);
            end
        end
        
        function plotUnitTgt(varargin)
            mesh = varargin{1};
            if (nargin == 2)
                spc = varargin{2};
            else
                spc = 'r';
            end
            switch mesh.type
                case 'segment'
                    Xctr = mesh.ctr;
                    Vnrm = mesh.tgt;
                    quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),...
                        'Color',spc);
                    
                otherwise
                    plotUnitTgt(mesh.edg,spc);
            end
        end
        
        function plotUnitNrm(varargin)
            mesh = varargin{1};
            assert(mesh.dim<=3);
            spc  = 'r';
            if (nargin == 2)
                spc = varargin{2};
            end
            Xctr = mesh.ctr;
            Vnrm = mesh.nrm;
            quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),...
                'Color',spc);
        end
        
        function plotNrmEdg(varargin)
            b = ~ishold;
            this = varargin{1};
            assert(this.dim<=3);
            spc  = 'r';
            J = 1:this.nelt;
            if (nargin >= 2)
                spc = varargin{2};
            end
            if (nargin ==3)
                J = varargin{3};
            end
            switch this.type
                case 'segment'
                    plotNrm(this);
                case 'triangle'
                    [A,B,C] = ABCD(this);
                    I{1} = (A(J,:)+B(J,:))/2; 
                    I{2} = (B(J,:)+C(J,:))/2; 
                    I{3} = (A(J,:)+C(J,:))/2;
                    nu = this.nrmEdg(J);
                    for d = 1:3
                        dp1 = mod(d,3)+1;
                        quiver3(I{dp1}(:,1),I{dp1}(:,2),I{dp1}(:,3),...
                            nu{d}(:,1),nu{d}(:,2),nu{d}(:,3),...
                        'color',spc,'AutoScaleFactor',.3);
                        hold on
                    end
                otherwise
                    error('unavailable case');  
            end
            if b
                % Release hold if it was in off mode before function call.
                hold off;
            end
        end
        
        function plotOn(mesh,func)
            % Plot a function func defined on the mesh
            mshPlotOn(mesh,func);
        end
        
        %% Access properties
        
        function[nel] = nelt(this)
            % Number of elements
            nel = size(this.elt,1);
        end
        
        % By convention, length(mesh) is the number of elements.
        function l = length(mesh)
            l = size(mesh.elt,1);
        end
        
        function[nv] = nvtx(this)
            % Nb of vtx
            nv = size(this.vtx,1);
        end
        
        function[nf] = nfce(this)
            assert(this.dim >=3)
            nf = length(this.fce);
        end
        
        function[ne] = nedg(this)
            assert(this.dim>=2)
            ne = length(this.edg);
        end
        
        function[nc] = ncol(this)
            % Number of colors
            nc = size(this.col,1);
        end
        
        function[d] = dim(this)
            % 4-> Tetrahedron mesh, 3-> Trimesh, 2-> edge mesh, 1-> particule
            % mesh.
            d = size(this.elt,2);
        end
        
        function[s] = type(this)
            switch this.dim
                case 1
                    s = 'point';
                case 2
                    s = 'segment';
                case 3
                    s = 'triangle';
                case 4
                    s = 'tetrahedron';
            end
        end
        
        % Local vertices of element. 
        function varargout = ABCD(this)
            % [A,B,C,D] = ABCD(mesh) : if mesh is a tetra mesh, A(i,:) contains the
            % coordinates of the first vtx of mesh.elts(i,:), B(i,:) the
            % second, and so on. For a particle, edge and triangle mesh, only
            % A, resp A,B, resp A,B,C, are defined.
            
            varargout = cell(1,this.dim);
            for i = 1:this.dim
                varargout{i} = this.vtx(this.elt(:,i),:);
            end
            for i = this.dim+1:4
                varargout{i} = NaN;
            end
        end
        
        % Size
        function s = size(varargin)
            s = size(varargin{1}.elt);
            if (nargin == 2)
                s = s(varargin{2});
            end
        end
        
        % Step
        function l = stp(mesh)
            mesh = mesh.edg;
            l    = mesh.vtx(mesh.elt(:,2),:) - mesh.vtx(mesh.elt(:,1),:);
            l    = sqrt(sum(l.^2,2));
            l    = [min(l) max(l) mean(l) std(l)];
        end
        
        % Number of space dimensions
        function b = is1d(mesh)
            b = is2d(mesh) && (max(abs(mesh.vtx(:,2))) < 1e-12);
        end
        
        function b = is2d(mesh)
            b = (max(abs(mesh.vtx(:,3))) < 1e-12);
        end
        
        % Mesh with no point / element
        function b = istrivial(mesh)
            b = mesh.length==0;
        end
        
        % Mesh without boundary.
        function b = isclosed(mesh)
            b = istrivial(mesh.bnd);
        end
        
        %% Set properties
        
        % Set color
        function mesh = color(mesh,c)
            assert(isnumeric(c));
            mesh.col(:) = c(:);
        end
        
        %% Geometric data
        
        % Center
        function X = ctr(this)
            X = zeros(length(this),3);
            [A{1},A{2},A{3},A{4}] = ABCD(this);
            for d = 1:this.dim
                X = X + (1/this.dim) .* A{d};
            end
        end
        
        % Nd-Volume
        function V = ndv(mesh)
            V = mshNdvolume(mesh);
        end
        
        % Unit normals
        function N = nrm(mesh)
            switch mesh.type
                case 'point'
                    error('No meaningful normal on point meshes');
                case 'segment'
                    if is2d(mesh)
                        % Take normal in the {Z = 0} plane
                        T = mesh.tgt;
                        N = T * [0 -1 0 ; 1 0 0 ; 0 0 0]';
                    else
                        error('Cannot define a normal on segment meshes in 3D');
                    end
                case 'triangle'
                    T = mesh.tgt;
                    N = cross(T{1},T{2},2);
                    N = N ./ (sqrt(sum(N.^2,2)) * [1 1 1]);
                case 'tetrahedron'
                    error('No meaningful normal on tetrahedron meshes');
            end
        end
        
        % Edge normals
        function Nu = nrmEdg(mesh,I)
            switch mesh.type
                case 'point'
                    error('No meaningful edge normal on point meshes')
                case 'segment'
                    Nu = mesh.nrm;
                    warning('Something harmful maybe here');
                    % Martin : I changed the definition.
                    % Previously, it was Nu{1} = mesh.tgt; Nu{2} = -mesh.tgt;
                case 'triangle'
                    Nu = cell(1,3);
                    t = mesh.tgt;
                    n = mesh.nrm;
                    for i = 1:3
                        ip1 = mod(i,3)+1;
                        Nu{i} = cross(t{ip1}(I,:),n(I,:),2);
                    end
                case 'tetrahedron'
                    error('No meaningful edge normal on tetrahedral meshes')
            end
        end
        
        % Tangents
        function T = tgt(this)
            switch this.type
                case 'point'
                    error('No meaningful tangent on point meshes')
                case 'segment'
                    [A,B] = ABCD(this);
                    T = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
                otherwise
                    T = cell(1,this.dim);
                    [A{1},A{2},A{3},A{4}] = ABCD(this);
                    for i = 1:this.dim
                        ip1 = mod(i,this.dim) + 1;
                        T{i} = (A{ip1} - A{i})./(sqrt(sum((A{ip1} - A{i}).^2,2))...
                            *[1 1 1]);
                        
                    end
            end
        end
        
        
        %% Submesh
        
        
        % Select elements
        function mesh = sub(mesh,Ielt)
            mesh.elt = mesh.elt(Ielt,:);
            mesh.col = mesh.col(Ielt);
            mesh     = clean(mesh);
        end
        
        % Select vertices geometrically 
        function mesh = subVtx(mesh,vtx)
            nvtx = size(vtx,1);
            mp1 = mesh.prt;
            mp2 = msh(vtx,(1:nvtx)');
            [~,I] = intersect(mp1,mp2);
            % I contains the indices of the vtx to be kept in mesh
            % Only keep elements that contain only the desired vertices:
            aux = sum(ismember(mesh.elt,I),2);
            J = find(aux==mesh.dim);
            mesh = sub(mesh,J);
        end
        
        % Face mesh
        function [mesh,elt2fce] = fce(mesh)
            assert(mesh.dim >= 3,'Segment and triangular meshes dont have faces');
            [mesh,elt2fce] = mshFace(mesh);
        end
        
        % Edges
        function [mesh,elt2edg] = edg(mesh)
            assert(mesh.dim >= 2,'Point meshes dont have edges');
            [mesh,elt2edg] = mshEdge(mesh);
        end
        
        % Particles
        function [mesh,elt2prt] = prt(mesh)
            [mesh,elt2prt] = mshParticle(mesh);
        end
        
        % Boundary
        function [mesh,I] = bnd(mesh)
            [mesh,I] = mshBoundary(mesh);
        end
        
        % Refine triangle mesh with midpoint
        function [mesh,Ir] = midpoint(varargin)
            mesh = varargin{1};
%             assert(strcmp(mesh.type,'triangle'));
            if (nargin == 1)
                I = (1:mesh.nelt)';
            else
                I = varargin{2};
            end
            switch mesh.type
                case 'segment'
                    [mesh,Ir] = mshMidpointSeg(mesh,I);
                case 'triangle'
                    [mesh,Ir] = mshMidpoint(mesh,I);
                otherwise
                    error('Unavailable case');
            end
        end
        
        % Perform several iterations of midpoint refinements on triangle
        % mesh
        function mesh = refine(mesh,spc)
            % The spc controls the refinement. It can be specified as
            % either
            % - A number : the refinement algorithm will stop well all
            % edges are smaller than spc
            % - An array of logicals of size nelt x 1 indicating which
            % elements are to be refined by midpoint
            % - An array of numbers speciying how many times each element
            % must be refined by the mipoint rule
            % - A function defined on R^3 such that f(X) is the number of
            % times an element lying at X should be refined.
            if ~exist('spc','var')||isempty(spc)
                spc = ones(this.nelt,1);
            end
            ref  = sum(mesh.ndv);
            if isnumeric(spc) && isequal(spc,1)
                mesh = midpoint(mesh);
            else
                mesh = mshRefine(mesh,spc);
            end
            % Check that volume is preserved.
            sol = sum(mesh.ndv);
            if (norm(ref-sol)/norm(ref) > 1e-12)
                error('msh.m : unavailable case')
            end
        end
        
        % HIERARCHICAL TREE
        function tree = tree(varargin)
            mesh = varargin{1};
            if (nargin == 1)
                typ = 'octree';
                Nlf = 1;
                fig = 0;
            elseif (nargin == 2)
                typ = varargin{2};
                Nlf = 1;
                fig = 0;
            elseif (nargin == 3)
                typ = varargin{2};
                Nlf = varargin{3};
                fig = 0;
            elseif (nargin == 4)
                typ = varargin{2};
                Nlf = varargin{3};
                fig = varargin{4};
            end
            tree = mshTree(mesh,typ,Nlf,fig);
        end
        
        
        %%  Algebra (Union, intersection, set difference)
        % Signature
        function M = sgn(mesh)
            % Creates a matrix, each row "describing" an element by giving
            % coordinates of its vertices and midpoint.
            M = cell(1,size(mesh.elt,2));
            for i = 1:length(M)
                M{i} = mesh.vtx(mesh.elt(:,i),:);
            end
            M = sort(cell2mat(M),2);
            M = [M,mesh.ctr];
            M = single(M);
        end
        
        % Remove duplicate vertices or elements.
        function [meshC,IA,IC] = unique(mesh)
            M         = sgn(mesh);
            [~,IA,IC] = unique(M,'rows','stable');
            meshC     = mesh.sub(IA);
        end
        
        % Intersection
        function [meshC,IA,IB] = intersect(meshA,meshB)
            assert(meshA.dim == meshB.dim,...
                'Intersection of meshes reserved for meshes of same type');
            A         = sgn(meshA);
            B         = sgn(meshB);
            [~,IA,IB] = intersect(A,B,'rows','stable');
            meshC     = meshA.sub(IA);
        end
        
        % DIFFERENCE
        function [meshC,IA] = setdiff(meshA,meshB)
            assert(meshA.dim == meshB.dim,...
                'Difference of meshes reserved for meshes of same type');
            A      = sgn(meshA);
            B      = sgn(meshB);
            [~,IA] = setdiff(A,B,'rows','stable');
            meshC  = meshA.sub(IA);
        end
        
        % Union
        function [meshC,IA,IB] = union(meshA,meshB)
            assert(meshA.dim == meshB.dim,...
                'Union of meshes reserved for meshes of same type');
            A         = sgn(meshA);
            B         = sgn(meshB);
            [~,IA,IB] = union(A,B,'rows','stable');
            meshA     = meshA.sub(IA);
            meshB     = meshB.sub(IB);
            vtxC      = [meshA.vtx ; meshB.vtx];
            eltC      = [meshA.elt ; meshB.elt + size(meshA.vtx,1)];
            colC      = [meshA.col ; meshB.col];
            meshC     = msh(vtxC,eltC,colC);
        end
        
        % Ismember
        function [IA,IB] = ismember(meshA,meshB)
            A       = sgn(meshA);
            B       = sgn(meshB);
            [IA,IB] = ismember(A,B,'rows');
        end
        
        % Image of mesh by a function
        function mesh = fct(mesh,fct)
            mesh.vtx = fct(mesh.vtx);
        end
        
        % Shuffle
        function mesh = shuffle(mesh,~)
            Nvtx = size(mesh.vtx,1);
            RPV = randperm(Nvtx);
            
            if nargin == 2
                warning('This was used');
            end
            % This seemed useless:
            % if (nargin == 2)
            %       RPV = (1:Nvtx)';
            % else
            %       RPV = randperm(Nvtx);
            % end
            mesh.vtx(RPV,:) = mesh.vtx;
            mesh.elt        = RPV(mesh.elt);
            RPE      = randperm(length(mesh));
            mesh.elt = mesh.elt(RPE,:);
            mesh.col = mesh.col(RPE,:);
        end
        
        % Equality (up to shuffle)
        function b = isequal(m1,m2,varargin)
            p = inputParser;
            p.addOptional('upToShuffle',true)
            p.parse(varargin{:});
            if ~p.Results.upToShuffle
                b = isequal(m1.vtx,m2.vtx) && isequal(m1.elt,m2.elt)...
                    && isequal(m1.col,m2.col);
            else
                if m1.dim ~= m2.dim
                    b = false;
                else
                    P1 = sgn(m1); P2 = sgn(m2);
                    P1 = sort(P1); P2 = sort(P2);
                    b = isequal(P1,P2);
                end
            end
            
        end
        
        
        %%  Geometric transformation
        
        % Translation
        function mesh = translate(mesh,U)
            mesh.vtx = mesh.vtx + ones(size(mesh.vtx,1),1)*U;
        end
        
        % Rotation
        function mesh = rotate(mesh,center,U,phi)
            % Rotate by angle phi around axis center + kU in trigonometric
            % direction. 
            N        = U./norm(U);
            mtemp = translate(mesh,-center);
            mtemp.vtx = cos(phi) * mtemp.vtx + ...
                (1-cos(phi)) .* ((mtemp.vtx*N')*N) + ...
                sin(phi) .* cross(ones(size(mtemp.vtx,1),1)*N,mtemp.vtx,2);
            mesh = translate(mtemp,center);
        end
        
        % Swap
        function mesh = swap(varargin)
            mesh = varargin{1};
            Ielt = 1:size(mesh.elt,1);
            if (nargin == 2)
                Ielt = varargin{2};
            end
            switch mesh.type
                case 'segment'
                    mesh.elt(Ielt,:) = mesh.elt(Ielt,[2 1]);
                case 'triangle'
                    mesh.elt(Ielt,:) = mesh.elt(Ielt,[2 1 3]);
            end
        end
        
        % Split
        function [mesh1,mesh2] = split(varargin)
            [mesh1,mesh2] = mshSplit(varargin{:});
        end
        
        % Explode
        function [new_m] = explode(m,rho)
            if ~exist('rho','var')||isempty(rho)
                rho = 0.92;
            end
            new_elt = zeros(length(m),dim(m));
            new_vtx = zeros(m.dim*m.nelt,3);
            nvtx = 0;
            
            for i = 1:length(m)
                vtx_i = m.vtx(m.elt(i,:),:);
                ctr = mean(vtx_i,1);
                delta_i = vtx_i - ones(m.dim,1)*ctr;
                
                new_vtx(m.dim*(i-1) + (1:m.dim),:) ...
                    = ctr + rho*delta_i;
                new_elt(i,:) = (nvtx+1):(nvtx+dim(m));
                nvtx= nvtx+dim(m);
            end
            
            new_m = msh(new_vtx,new_elt,m.col);
        end
        
    end
end
