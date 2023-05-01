% Post processing the Neumann Data obtained from solve_superconductor which
% lies in RWG space

function out = plot_field_magnet(TdA,TnA,bndmesh,J,omega_src,mu0,interior)
    Gamma = dom(bndmesh,3);
    plot(bndmesh);
    hold on;
    plot(omega_src.msh);

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

    %scatter3(X,Y,Z);
    %hold on;

    eval_pts = [X,Y,Z];
    interior_pts = interior(eval_pts);

    eval_pts = eval_pts(interior_pts == 0,:);

    Bfield = curl_DLP_Vectorial(Gamma,TdA,eval_pts)...
        -curl_SLP_Vectorial(Gamma,TnA,eval_pts)...
        + mu0* compute_vecpot_curl(J,omega_src,eval_pts);

    quiver3wrapper(eval_pts, Bfield, 'red');

    xlabel('x');
    ylabel('y');
    zlabel('z');

    out =[];
%     maxSrc = max(omega_src.msh.vtx);
%     minSrc = min(omega_src.msh.vtx);
%     meanSrc = 0.5(maxSrc+minSrc);

    
end