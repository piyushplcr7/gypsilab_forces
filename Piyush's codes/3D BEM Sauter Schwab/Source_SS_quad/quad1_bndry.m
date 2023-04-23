function [X, W] = quad1_bndry(Xqud, Wqud)


        n1 = Xqud(:, 1);
        n2 = Xqud(:, 2);
        n3 = Xqud(:, 3);
        z  = Xqud(:, 4);

       
        oneColumn = 0*z+1;
        
        Jac = Wqud .* z.^3 .* n1.^2 .* n2;
        
        W1 = Jac;
        X1 = z .* [oneColumn,...
            1-n1+n1.*n2,...
            1-n1.*n2.*n3,...
            1-n1];
        
        W2 = Jac;
        X2 = z .* [1-n1.*n2.*n3,...
            1-n1,...
            oneColumn,...
            1-n1+n1.*n2];
        
        W3 = Jac;
        X3 = z .* [oneColumn,...
            n1.*(1-n2+n2.*n3),...
            1-n1.*n2,...
            n1.*(1-n2)];
        
        W4 = Jac;
        X4 = z .* [1-n1.*n2,...
            n1.*(1-n2),...
            oneColumn,...
            n1.*(1-n2+n2.*n3)];
        
        W5 = Jac;
        X5 = z .* [1-n1.*n2.*n3,...
            n1.*(1-n2.*n3),...
            oneColumn,...
            n1.*(1-n2)];
        
        W6 = Jac;
        X6 = z .* [oneColumn,...
            n1.*(1-n2),...
            1-n1.*n2.*n3,...
            n1.*(1-n2.*n3)];
        
        W = [W1;W2;W3;W4;W5;W6];

        X = [X1;X2;X3;X4;X5;X6];

end