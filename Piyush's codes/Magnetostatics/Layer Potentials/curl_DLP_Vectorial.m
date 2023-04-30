function out = curl_DLP_Vectorial(Gamma,TdA,eval_pts)

    RWG = fem(Gamma.msh,'RWG');
    % using nxNED = RWG to evaluate Rgy at quadrature points
    Rgy = reconstruct(TdA,Gamma,RWG);

    % Getting the quadrature weights and nodes
    [Y,W] = Gamma.qud;
    NY = size(Y,1);
    NX = size(eval_pts,1);

    % creating pairs of all quadrature points and eval_pts
    XX = repelem(eval_pts,NY,1);
    YY = repmat(Y,NX,1);
    RgYY = repmat(Rgy,NX,1);

    % Kernel gradx gradx G
    integrand = 3/(4*pi) * (XX-YY).*dot(XX-YY,RgYY,2)./vecnorm(XX-YY,2,2).^5 ...
                - 1/(4*pi) * RgYY./vecnorm(XX-YY,2,2).^3;
    
    integrands = reshape(integrand,[NY 3*NX]);

    integrals = sum(W.*integrands,1);

    out = reshape(integrals,[NX 3]);
end