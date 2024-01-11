function [potX,gradPotX] = computeNewtonPotentialAndGrad(X,omega,rho)
    % Potential to be computed at X
    
    % Quadrature points over integration domain
    [Y,W] = omega.qud;

    % Size of point list
    NX = size(X,1);
    NY = size(Y,1);

    % Creating combinations of X and Y. repmat for Y, repelem for X
    XX = repelem(X,NY,1);
    YY = repmat(Y,NX,1);

    % Fundamental solution and its gradient
    Gxy = 1/4/pi./vecnorm(XX-YY,2,2);
    gradxGxy = 1/4/pi * (YY-XX)./vecnorm(XX-YY,2,2).^3;

    % Reshaping to get G in the form
    % [G(X1,Y) G(X2,Y) ....]
    Gxy = reshape(Gxy,NY,NX);

    % Reshaping gradxG 
    gradxGxy = reshape(gradxGxy,NY,3*NX);

    % Computing the integrals
    rhoY = rho(Y);
    potX = sum(W.*Gxy.*rhoY,1)';
    gradPotX = sum(W.*gradxGxy.*rhoY,1);
    gradPotX = reshape(gradPotX,NX,3);

end