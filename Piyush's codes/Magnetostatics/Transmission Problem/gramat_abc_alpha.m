% Gramian matrix calculator for the order in which fields were used
function gramat = gramat_abc_alpha()
    
    l = 6*pi;
    N = 5000;
    mesh = mshCube(N,[l l l]);
    
    omega = dom(mesh,4);
    [X,W] = omega.qud;
    
    % indices for velocity fields go from 0 to kappa-1
    kappa = 3;
    
    % Map of velocity field to index for x pointing fields
    abc_alpha = zeros(kappa^3,3);
    
    for a = 0:kappa-1
        for b = 0:kappa-1
            for c = 0:kappa-1
                idx = a + kappa * b + kappa^2 * c + 1;
                abc_alpha(idx,:) = [a b c];
            end
        end
    end
    
    gramat_small = zeros(kappa^3,kappa^3);
    gramat = zeros(3*kappa^3,3*kappa^3);
    manually_constructed = gramat_small;
    
    for i = 1:kappa^3
        a1 = abc_alpha(i,1);
        b1 = abc_alpha(i,2);
        c1 = abc_alpha(i,3);
        [Vel1, DVel1] = getCosVelDVel(a1,b1,c1,1);
    
        for j = 1:kappa^3
            a2 = abc_alpha(j,1);
            b2 = abc_alpha(j,2);
            c2 = abc_alpha(j,3);
            [Vel2, DVel2] = getCosVelDVel(a2,b2,c2,1);
    
%             % L2 inner pdt
%             l2pdt = sum(W.*dot(Vel1(X),Vel2(X),2),1);
%     
%             % H1 inner pdt
%             h1pdt = sum(W.*(dot(DVel1{1}(X),DVel2{1}(X),2)...
%                            +dot(DVel1{2}(X),DVel2{2}(X),2)...
%                            +dot(DVel1{3}(X),DVel2{3}(X),2)),1);
%     
%             gramat_small(i,j) = l2pdt + h1pdt;
            
            % Manually constructed l2 inner product
            manually_constructed(i,j) = Ia1a2(a1,a2,l/2/pi)...
                                        * Ia1a2(b1,b2,l/2/pi)...
                                        * Ia1a2(c1,c2,l/2/pi);
    
            % Manually constructed h1 inner product
            manually_constructed(i,j) = manually_constructed(i,j) + a1*a2*Ja1a2(a1,a2,l/2/pi)*Ia1a2(b1,b2,l/2/pi)*Ia1a2(c1,c2,l/2/pi)...
                                       +b1*b2*Ia1a2(a1,a2,l/2/pi)*Ja1a2(b1,b2,l/2/pi)*Ia1a2(c1,c2,l/2/pi)...
                                       +c1*c2*Ia1a2(a1,a2,l/2/pi)*Ia1a2(b1,b2,l/2/pi)*Ja1a2(c1,c2,l/2/pi);
    
        end
    
    end
    
%     gramat(1:kappa^3,1:kappa^3) = gramat_small;
%     gramat(kappa^3+1:2*kappa^3,kappa^3+1:2*n^3) = gramat_small;
%     gramat(2*kappa^3+1:3*kappa^3,2*kappa^3+1:3*kappa^3) = gramat_small;

    gramat(1:kappa^3,1:kappa^3) = manually_constructed;
    gramat(kappa^3+1:2*kappa^3,kappa^3+1:2*kappa^3) = manually_constructed;
    gramat(2*kappa^3+1:3*kappa^3,2*kappa^3+1:3*kappa^3) = manually_constructed;

end