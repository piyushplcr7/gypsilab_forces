function [X, W] = quad3_bndry(Xqud, Wqud)


        n1 = Xqud(:, 1);
        n2 = Xqud(:, 2);
        n3 = Xqud(:, 3);
        z  = Xqud(:, 4);
        
        oneColumn = 0*z+1;
        
        
        Jac1 = Wqud .* z.^3 .* n2;
        
        W1 = Jac1;
        X1 = z .* [oneColumn,...
            n1,...
            n2,...
            n2.*n3];
        
        W2 = Jac1;
        X2 = z .* [n2,...
            n2.*n3,...
            oneColumn,...
            n1];
        
        
        W = [W1;W2];

        X = [X1;X2];

end