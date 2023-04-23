function [X, W] = quad2_bndry(Xqud, Wqud)


        n1 = Xqud(:, 1);
        n2 = Xqud(:, 2);
        n3 = Xqud(:, 3);
        z  = Xqud(:, 4);
        
        oneColumn = 0*z+1;
        
        
        Jac1 = Wqud .* z.^3 .* n1.^2 .* n2;
        
        W1 = Jac1;
        X1 = z .* [oneColumn,...
            n1,...
            1-n1.*n2.*n3,...
            n1.*n2.*(1-n3)];
        
        W2 = Jac1;
        X2 = z .* [1-n1.*n2,...
            n1.*(1-n2),...
            oneColumn,...
            n1.*n2.*n3];
        
        W3 = Jac1;
        X3 = z .* [1-n1.*n2.*n3,...
            n1.*n2.*(1-n3),...
            oneColumn,...
            n1];
        
        W4 = Jac1;
        X4 = z .* [1-n1.*n2.*n3,...
            n1.*(1-n2.*n3),...
            oneColumn,...
            n1.*n2];
        
        
        Jac2 = Wqud .* z.^3 .* n1.^2;
        W5 = Jac2;
        X5 = z .* [oneColumn,...
            n1.*n3,...
            1-n1.*n2,...
            n1.*(1-n2)];
        
        W = [W1;W2;W3;W4;W5];

        X = [X1;X2;X3;X4;X5];

end