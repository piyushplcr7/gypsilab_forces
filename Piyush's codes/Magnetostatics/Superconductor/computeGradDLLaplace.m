function out = computeGradDLLaplace(bndmesh,Tdu,eval_pts)
    % BEM Space for Tdu
    P1 = fem(bndmesh,'P1');
    Gamma = dom(bndmesh,3);

    % Computing psi at quadrature points
    g = reconstruct(Tdu,Gamma,P1);

    N_eval_pts = size(eval_pts,1);
    
    [Y,W] = Gamma.qud;
    normals = Gamma.qudNrm;
    ng = g.*normals;

    N_qud_pts = size(Y,1);

    eval_pts_ex = repelem(eval_pts,N_qud_pts,1); % evaluation points X
    Y_ex = repmat(Y,N_eval_pts,1); % integration variable Y
    % X-Y
    xmy = eval_pts_ex - Y_ex;

    ng = repmat(ng,N_eval_pts,1);

    xmydotng = dot(xmy,ng,2);

    vectorintegrand = -3/(4*pi).*xmy.*xmydotng./vecnorm(xmy,2,2).^5 + 1/(4*pi)*ng./vecnorm(xmy,2,2).^3;
    vectorintegrand1 = reshape(vectorintegrand(:,1),[N_qud_pts,N_eval_pts]);
    vectorintegrand2 = reshape(vectorintegrand(:,2),[N_qud_pts,N_eval_pts]);
    vectorintegrand3 = reshape(vectorintegrand(:,3),[N_qud_pts,N_eval_pts]);

    t1 = sum(W.*g.*vectorintegrand1)';
    t2 = sum(W.*g.*vectorintegrand2)';
    t3 = sum(W.*g.*vectorintegrand3)';

    out = [t1 t2 t3];

end