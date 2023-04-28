function [X, W] = quad2_cross(Xqud, Wqud)

        n1 = Xqud(:, 1);
        n2 = Xqud(:, 2);
        n3 = Xqud(:, 3);
        z1 = Xqud(:, 4);
        z2 = Xqud(:, 5);
        
        
        oneColumn = 0*z1+1;
        
        Jac1 = Wqud .* z1.^4 .* z2.^3 .*n2; 
        
        W1 = Jac1;
        X1 = (1-z1) .* [oneColumn, 0*oneColumn, oneColumn, 0*oneColumn, 0*oneColumn] + ... 
                 z1 .* z2.* [n2.*n3,...
                             n2.*(1-n3),...
                             0*oneColumn,...
                             n1,...
                             (1-n1)];

        W2 = Jac1;
        X2 = (1-z1) .* [oneColumn, 0*oneColumn, oneColumn, 0*oneColumn, 0*oneColumn] + ... 
                 z1 .* z2.* [n1,...
                             (1-n1),...
                             0*oneColumn,...
                             n2.*n3,...
                             n2.*(1-n3)];
        
        W3 = Jac1;                 
        X3 = (1-z1) .* [oneColumn, 0*oneColumn, oneColumn, 0*oneColumn, 0*oneColumn] + ... 
                       z1 .* z2.* [0*oneColumn,...
                                   n1,...
                                   n2.*n3,...
                                   n2.*(1-n3),...
                                   (1-n2)];
        
        
                               
        Jac2 = Wqud .* z1.^4 .* z2.^3 .* n1.^2 .* n2;
        
        W4 = Jac2;
        X4 = (1-z1) .* [oneColumn, 0*oneColumn, oneColumn, 0*oneColumn, 0*oneColumn] + ... 
                       z1 .* z2.* [0*oneColumn,...
                                   oneColumn,...
                                   n1.*n2.*n3,...
                                   n1.*n2.*(1-n3),...
                                   n1.*(1-n2)];
                               
       W = [W1;W2;W3;W4];
       X = [X1;X2;X3;X4];
        

end