classdef msFem < fem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function[] = surf(this,val)
            [~,~,restrictions] = dofMS(this);
            typ = this.typ(7:end);
            for i = 1:this.msh.npanels
                pani = this.msh.panels{i};
                Nrm = pani.nrm;
                n = mean(Nrm,1);
                pani = translate(pani,-1e-3*n);
                femi = fem(pani,typ);
                surf(femi,restrictions{i}*val);
                hold on;
            end
        end
    end
end

