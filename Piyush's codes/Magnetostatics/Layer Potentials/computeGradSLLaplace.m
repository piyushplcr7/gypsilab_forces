function out = computeGradSLLaplace(bndmesh,Tnu,eval_pts)
    % BEM Space for Tnu
    P0 = fem(bndmesh,'P0');
    Gamma = dom(bndmesh,3);

    % Computing psi at quadrature points
    psi = reconstruct(Tnu,Gamma,P0);

    N_eval_pts = size(eval_pts,1);
    
    [Y,W] = Gamma.qud;
    N_qud_pts = size(Y,1);

    eval_pts_ex = repelem(eval_pts,N_qud_pts,1); % evaluation points X
    Y_ex = repmat(Y,N_eval_pts,1); % integration variable Y
    % Y - X
    ymx = Y_ex - eval_pts_ex; 

    gradxG = 1/(4*pi) * ymx ./vecnorm(ymx,2,2).^3; %(Matrix with 3 columns)

    gradxGmat = cell(3,1);
    gradxGmat{1} = reshape(gradxG(:,1),[N_qud_pts,N_eval_pts]);
    gradxGmat{2} = reshape(gradxG(:,2),[N_qud_pts,N_eval_pts]);
    gradxGmat{3} = reshape(gradxG(:,3),[N_qud_pts,N_eval_pts]);

    t1 = sum(W.*psi.*gradxGmat{1})';
    t2 = sum(W.*psi.*gradxGmat{2})';
    t3 = sum(W.*psi.*gradxGmat{3})';

    out = [t1 t2 t3];

end