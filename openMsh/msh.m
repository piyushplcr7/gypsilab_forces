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
            this = clean(this);
        end
        
        
        % CLEAN % Remove duplicate vertices, remove unused vertices and relabel
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
                    Xctr = mesh.ctr;
                    Vtgt = mesh.tgt;
                    Vtgt = Vtgt.*(mesh.ndv*([1 1 1]));
                    quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vtgt(:,1),Vtgt(:,2),Vtgt(:,3),...
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
            assert(mesh.dim <=2);
            spc  = 'r';
            if (nargin == 2)
                spc = varargin{2};
            end
            Xctr = mesh.ctr;
            Vnrm = mesh.tgt;
            quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),spc);
            
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
            quiver3(Xctr(:,1),Xctr(:,2),Xctr(:,3),Vnrm(:,1),Vnrm(:,2),Vnrm(:,3),spc);
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
        
        function[ne] = nedg(this)
            assert(this.dim>=2)
            ne = length(this.edg);
        end
        
        function[nf] = nfce(this)
            assert(this.dim >=3)
            nf = length(this.fce);
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
        
        function varargout = ABCD(this)
            % [A,B,C,D] = ABCD(mesh) : if mesh is a tetra mesh, A(i,:) contains the
            % coordinates of the first vtx of mesh.elts(i,:), B(i,:) the
            % second, and so on. For a particle, edge and triangle mesh, only
            % A, resp A,B, resp A,B,C, are defined.
            
            if nargout > this.dim
                error('too many output arguments');
            end
            varargout = cell(1,this.dim);
            for i = 1:this.dim
                varargout{i} = this.vtx(this.elt(:,i),:);
            end
        end
        
        % SIZE
        function s = size(varargin)
            s = size(varargin{1}.elt);
            if (nargin == 2)
                s = s(varargin{2});
            end
        end
        
        % STEP
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
        
        %% Geometric data
        
        % CENTER
        function X = ctr(mesh)
            X = zeros(size(mesh.elt,1),size(mesh.vtx,2));
            for i = 1:size(mesh.elt,2)
                X = X + (1/size(mesh.elt,2)) .* mesh.vtx(mesh.elt(:,i),:);
            end
        end
        
        % ND-VOLUME
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
        function Nu = nrmEdg(mesh)
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
                    for i = 1:3
                        Nu{i} = cross(mesh.tgt{i},mesh.nrm,2);
                    end
                case 'tetrahedron'
                    error('No meaningful edge normal on tetrahedral meshes')
            end
        end
        
        % Tangents
        function T = tgt(mesh)
            switch mesh.type
                case 'point'
                    error('No meaningful tangent on point meshes')
                case 'segment'
                    [A,B] = ABCD(mesh);
                    T = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
                case 'triangle'
                    % Return a cell containing {AB,BC,CA} where A,B,C is the
                    % order of the vertices as they appear in mesh.elt
                    T = cell(1,3);
                    [A,B,C] = ABCD(mesh);
                    T{1} = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
                    T{2} = (C-B)./(sqrt(sum((C-B).^2,2))*[1 1 1]);
                    T{3} = (A-C)./(sqrt(sum((A-C).^2,2))*[1 1 1]);
                    % for i = 1:3
                    %    ip1  = mod(i,3)+1;
                    %    ip2  = mod(ip1,3)+1;
                    %    A    = mesh.vtx(mesh.elt(:,ip1),:);
                    %    B    = mesh.vtx(mesh.elt(:,ip2),:);
                    %    T{i} = (B-A)./(sqrt(sum((B-A).^2,2))*[1 1 1]);
            end
        end
        
        
        % SWAP
        function mesh = swap(varargin)
            mesh = varargin{1};
            Ielt = 1:size(mesh.elt,1);
            if (nargin == 2)
                Ielt = varargin{2};
            end
            if (size(mesh,2) == 2)
                mesh.elt(Ielt,:) = mesh.elt(Ielt,[2 1]);
            elseif (size(mesh,2) == 3)
                mesh.elt(Ielt,:) = mesh.elt(Ielt,[2 1 3]);
            else
                error('msh.m : unavailable case')
            end
        end
        
        % COLOURS
        function mesh = color(mesh,c)
            mesh.col(:) = c;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBMESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SUBMESHING
        function mesh = sub(mesh,Ielt)
            mesh.elt = mesh.elt(Ielt,:);
            mesh.col = mesh.col(Ielt);
            mesh     = clean(mesh);
        end
        
        % FACES
        function [mesh,elt2fce] = fce(mesh)
            [mesh,elt2fce] = mshFace(mesh);
        end
        
        % EDGES
        function [mesh,elt2edg] = edg(mesh)
            [mesh,elt2edg] = mshEdge(mesh);
        end
        
        % PARTICLES
        function [mesh,elt2prt] = prt(mesh)
            [mesh,elt2prt] = mshParticle(mesh);
        end
        
        % BOUNDARY
        function mesh = bnd(mesh)
            mesh = mshBoundary(mesh);
        end
        
        % MIDPOINT
        function [mesh,Ir] = midpoint(varargin)
            mesh = varargin{1};
            if (nargin == 1)
                I = (1:mesh.nelt)';
            else
                I = varargin{2};
            end
            [mesh,Ir] = mshMidpoint(mesh,I);
        end
        
        % REFINE WITH MIDPOINT
        function mesh = refine(varargin)
            mesh = varargin{1};
            ref  = sum(mesh.ndv);
            if (nargin == 1)
                mesh = midpoint(mesh);
            else
                mesh = mshRefine(mesh,varargin{2});
            end
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ALGEBRA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SIGNATURE
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
        
        % UNICITY
        function [meshC,IA,IC] = unique(mesh)
            M         = sgn(mesh);
            [~,IA,IC] = unique(M,'rows','stable');
            meshC     = mesh.sub(IA);
        end
        
        % INTERSECTION
        % Non-symmetric operation
        function [meshC,IA,IB] = intersect(meshA,meshB)
            A         = sgn(meshA);
            B         = sgn(meshB);
            [~,IA,IB] = intersect(A,B,'rows','stable');
            meshC     = meshA.sub(IA);
        end
        
        % DIFFERENCE
        function [meshC,IA] = setdiff(meshA,meshB)
            A      = sgn(meshA);
            B      = sgn(meshB);
            [~,IA] = setdiff(A,B,'rows','stable');
            meshC  = meshA.sub(IA);
        end
        
        % UNION
        function [meshC,IA,IB] = union(meshA,meshB)
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
        
        % ISMEMBER
        function [IA,IB] = ismember(meshA,meshB)
            A       = sgn(meshA);
            B       = sgn(meshB);
            [IA,IB] = ismember(A,B,'rows');
        end
        
        % FUNCTION
        function mesh = fct(mesh,fct)
            mesh.vtx = fct(mesh.vtx);
        end
        
        % SHUFFLE
        function mesh = shuffle(mesh,varargin)
            Nvtx = size(mesh.vtx,1);
            if (nargin == 2)
                RPV = (1:Nvtx)';
            else
                RPV = randperm(Nvtx);
            end
            mesh.vtx(RPV,:) = mesh.vtx;
            mesh.elt        = RPV(mesh.elt);
            RPE      = randperm(length(mesh));
            mesh.elt = mesh.elt(RPE,:);
            mesh.col = mesh.col(RPE,:);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% TRANSFORMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TRANSLATION
        function mesh = translate(mesh,U)
            mesh.vtx = mesh.vtx + ones(size(mesh.vtx,1),1)*U;
        end
        
        % ROTATION
        function mesh = rotate(mesh,U,phi)
            N        = U./norm(U);
            mesh.vtx = cos(phi) * mesh.vtx + ...
                (1-cos(phi)) .* ((mesh.vtx*N')*N) + ...
                sin(phi) .* cross(ones(size(mesh.vtx,1),1)*N,mesh.vtx,2);
        end
        
        % SPLIT
        function [mesh1,mesh2] = split(varargin)
            [mesh1,mesh2] = mshSplit(varargin{:});
        end
        
        function [new_m] = explode(m,rho)
            if ~exist('rho','var')||isempty(rho)
                rho = 0.92;
            end
            new_elt = zeros(length(m),dim(m));
            new_vtx = [];
            nvtx = 0;
            
            for i = 1:length(m)
                vtx_i = m.vtx(m.elt(i,:),:);
                ctr = mean(vtx_i,1);
                delta_i = vtx_i - ones(m.dim,1)*ctr;
                
                new_vtx = [new_vtx; ctr + rho*delta_i];
                new_elt(i,:) = (nvtx+1):(nvtx+dim(m));
                nvtx= nvtx+dim(m);
            end
            
            new_m = msh(new_vtx,new_elt,m.col);
        end
        
    end
end
