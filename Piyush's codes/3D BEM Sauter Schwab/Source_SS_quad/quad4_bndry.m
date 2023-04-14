
function [X, W] = quad4_bndry(Xqud, Wqud)


        n1 = Xqud(:, 1);
        n2 = Xqud(:, 2);
        n3 = Xqud(:, 3);
        n4 = Xqud(:, 4);
        
        X = [n1,...
            n2 .* n1,...
            n3,...
            n4 .* n3];
        
        W = Wqud .* n1 .* n3;
        
end