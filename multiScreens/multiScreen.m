classdef multiScreen
    
    properties
        panels; % List of meshes of each Gamma_j
        vtx;
        I
    end
    
    methods
        
        % Class Constructor
        
        function[this] = multiScreen(pan)
            if nargin == 0
                this.panels = cell(0);
            else
                 this.panels = pan;
                 vtx = [];
                 for i = 1:length(pan)
                     vtx = [vtx; pan{i}.vtx];
                 end
                 [this.vtx,this.I] = unique(vtx,'rows');
            end
        end
        
        % Methods
        function np = npanels(this)
            np = length(this.panels);
        end
        
        
        function[] = plot(varargin)
            this = varargin{1};
            plot(this.panels{1},varargin{2:end});
            hold on
            for i = 2:this.npanels
                plot(this.panels{i},varargin{2:end})
            end
        end
        
        function[this] = rotate(this,center,U,phi)
            for i = 1:this.npanels
                this.panels{i} = rotate(this.panels{i},center,U,phi);
            end
        end
        function[E] = elt(this)
            E = this.singleMesh.elt;
        end
        function[sg] = singleMesh(this)
            sg = this.panels{1};
            for i = 2:this.npanels
                sg = union(sg,this.panels{i});
            end
            sg = mshClean(sg);
        end
        function[b] = bnd(this)
            b = bnd(this.singleMesh);
        end
        function[M] = regularize(varargin)
            
        end
        
        
    end
end

