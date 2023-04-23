% curl SLP vectorial

function curlA = curl_SLP_Vectorial(omega_src,psi_coeffs,X)
     % Green's kernel
    GradGy = cell(3,1);
    GradGy{1} = @(X,Y)1/4/pi*femGreenKernel(X,Y,'grady[1/r]1',0);
    GradGy{2} = @(X,Y)1/4/pi*femGreenKernel(X,Y,'grady[1/r]2',0);
    GradGy{3} = @(X,Y)1/4/pi*femGreenKernel(X,Y,'grady[1/r]3',0);

    % Getting the quadrature weights and nodes
    [X_src,W_src,~] = omega_src.qud;

    % Evaluating the source current at source quadrature points
    % Evaluating source density at source quadrature nodes
    % BEM space in which the source density lies
    RWG = fem(omega_src.msh,'RWG');
    J_X_src= reconstruct(psi_coeffs,omega_src,RWG);
    N_qud_pts = size(X,1);
    N_qud_pts_src = size(X_src,1);
    
    % Repeating the source quadrature points N_qud_pts times
    X_src_extended = repmat(X_src,N_qud_pts,1);
    X_extended = repelem(X,N_qud_pts_src,1);   
    
    grad_green_fn_1 = GradGy{1}(X_extended,X_src_extended);
    grad_green_fn_2 = GradGy{2}(X_extended,X_src_extended);
    grad_green_fn_3 = GradGy{3}(X_extended,X_src_extended);

    % Reshaping
    grad_green_fn_1 = reshape(grad_green_fn_1,[N_qud_pts_src,N_qud_pts]);
    grad_green_fn_2 = reshape(grad_green_fn_2,[N_qud_pts_src,N_qud_pts]);
    grad_green_fn_3 = reshape(grad_green_fn_3,[N_qud_pts_src,N_qud_pts]);

    % Computing the cross product in the integrand
    blabla1 = J_X_src(:,2).*grad_green_fn_3 - J_X_src(:,3).*grad_green_fn_2;
    blabla2 = J_X_src(:,3).*grad_green_fn_1 - J_X_src(:,1).*grad_green_fn_3;
    blabla3 = J_X_src(:,1).*grad_green_fn_2 - J_X_src(:,2).*grad_green_fn_1;

    % Computing the components of the integral, in the end it contains the curl
    % of the vector potential at the quadrature nodes of Gamma
    blabla1 = sum(W_src.*blabla1,1)';
    blabla2 = sum(W_src.*blabla2,1)';
    blabla3 = sum(W_src.*blabla3,1)';
    
    curlA = [blabla1 blabla2 blabla3];


end