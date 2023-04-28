function [X, W] = quad4_cross(Xqud, Wqud)


        n1 = Xqud(:, 1);
        n2 = Xqud(:, 2);
        n3 = Xqud(:, 3);
        n4 = Xqud(:, 4);
        n5 = Xqud(:, 5);
        
        


        W1 = Wqud .* n1.^2 .* n2 .* n4;
        X1 =[n4.*n5,...
             n4.*(1-n5),...
             n1.*n2.*n3,...
             n1.*n2.*(1-n3),...
             n1.*(1-n2)];
        
        W = W1;
        X = X1;


end