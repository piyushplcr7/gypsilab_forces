% Post processing the Neumann Data obtained from solve_superconductor which
% lies in RWG space

function out = plot_field(TnA,bndmesh,J,omega_src)
    Gamma = dom(bndmesh,3);
    plot(bndmesh);
    hold on;
    %plot(omega_src.msh);

    maxGamma = max(bndmesh.vtx);
    minGamma = min(bndmesh.vtx);
    meanGamma = 0.5*(maxGamma+minGamma);
    newMinGamma = 2 * minGamma - meanGamma;
    newMaxGamma = 2 * maxGamma - meanGamma;

    N = 20;
    x = linspace(newMinGamma(1),newMaxGamma(1),N);
    y = linspace(newMinGamma(2),newMaxGamma(2),N);

    [X,Y] = meshgrid(x,y);
    X = reshape(X,[N*N,1]);
    Y = reshape(Y,[N*N,1]);
    Z = meanGamma(3) * ones(N*N,1);

    eval_pts = [X,Y,Z];

    Bfield = -curl_SLP_Vectorial(Gamma,TnA,eval_pts) + compute_vecpot_curl(J,omega_src,eval_pts);

    quiver3wrapper(eval_pts, Bfield, 'red');

    out =[];
%     maxSrc = max(omega_src.msh.vtx);
%     minSrc = min(omega_src.msh.vtx);
%     meanSrc = 0.5(maxSrc+minSrc);

    
end