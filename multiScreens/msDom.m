classdef msDom < dom
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Class constructor
        
        function this = msDom(varargin)
            if nargin == 0
            else
                this.msh = varargin{1};
                this.gss = varargin{2};
            end
        end
        
        
        function [X,W,I] = qud(this)
            X = []; W = [];
            gs = this.gss;
            ms = this.msh;
            I = cell(ms.npanels,1);
            start = 0;
            for i = 1:ms.npanels
                pani = ms.panels{i};
                domi = dom(pani,gs);
                [x,w] = qud(domi);
                X = [X;x]; W = [W;w];
                I{i} = (start+1):(start+length(w));
                start = start+length(w);
            end
        end
    end
end

