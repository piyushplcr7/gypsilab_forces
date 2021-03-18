classdef Florian
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        adress
        tel
    end
    
    methods
        function obj = getAdress(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = getTel(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


% S = singleLayer(m);
