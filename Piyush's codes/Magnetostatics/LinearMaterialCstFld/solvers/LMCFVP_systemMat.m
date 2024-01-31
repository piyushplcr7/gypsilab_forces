% Solver linear material constant field vector potential
% Solve Transmission Problem. Returns the exterior traces 

function blockopr = LMCFVP_systemMat(bndmesh_i,bndmesh_e,mu,mu_0)
    %% BEM Spaces
    Gamma_i = dom(bndmesh_i,3);

    Gamma_e = dom(bndmesh_e,3);

    % BEM spaces
    % curl conforming -> Dirichlet trace
    NED_i = fem(bndmesh_i,'NED'); 
    P1_i = fem(bndmesh_i,'P1');
    % Div conforming with div0 constraint -> Neumann trace
    DIV0_i = nxgrad(P1_i); 
    % Kernel of the surface curl operator
    Ker_curl_i = grad(P1_i); 
    % Div conforming space 
    DIV_i = fem(bndmesh_i,'RWG');

    P1_e = fem(bndmesh_e,'P1');
    DIV0_e = nxgrad(P1_e);

    %% Galerkin matrices 
    % For operator A
    Aii = single_layer(Gamma_i,DIV0_i,DIV0_i);
    % For operator C
    Cii = double_layer_magnetostatics(Gamma_i,DIV0_i,DIV_i);
    % For operator B, we can simply use C'
    Bii = Cii';
    % For operator N
    Nii = -single_layer(Gamma_i,DIV_i.div,DIV_i.div);

    % vector to enforce zero mean for P1 functions
    veci = integral(Gamma_i,P1_i);

    ortho = single_layer(Gamma_i,NED_i,Ker_curl_i); % Uses my implementation.

    % Cross matrices
    Aei = single_layer_cross(Gamma_i,Gamma_e,DIV0_i,DIV0_e);
    Cie = double_layer_magnetostatics_cross(Gamma_e,Gamma_i,DIV0_e,DIV_i);

    Aee = single_layer(Gamma_e,DIV0_e,DIV0_e);
    vece = integral(Gamma_e,P1_e);

    %% Linear System

    %                 P1_i.ndof,          DIV_i.ndof,                  P1_e.ndof,                   1,                  P1_i.ndof,                   1,                   1; 
    blockopr = [(1+mu/mu_0)*Aii,               2*Cii,                        Aei,                veci, zeros(P1_i.ndof,P1_i.ndof),  zeros(P1_i.ndof,1),  zeros(P1_i.ndof,1); %P1_i.ndof
                         -2*Bii,    -(1+mu_0/mu)*Nii,                       -Cie', zeros(DIV_i.ndof,1),                      ortho, zeros(DIV_i.ndof,1), zeros(DIV_i.ndof,1); % DIV_i.ndof
                           Aei',                 Cie,                        Aee,  zeros(P1_e.ndof,1), zeros(P1_e.ndof,P1_i.ndof),                vece,  zeros(P1_e.ndof,1); % P1_e.ndof
                          veci', zeros(1,DIV_i.ndof),         zeros(1,P1_e.ndof),                   0,         zeros(1,P1_i.ndof),                   0,                   0; % 1
     zeros(P1_i.ndof,P1_i.ndof),              ortho', zeros(P1_i.ndof,P1_e.ndof),  zeros(P1_i.ndof,1), zeros(P1_i.ndof,P1_i.ndof),  zeros(P1_i.ndof,1),                veci; % P1_i.ndof
             zeros(1,P1_i.ndof), zeros(1,DIV_i.ndof),                      vece',                   0,         zeros(1,P1_i.ndof),                   0,                   0; %1
             zeros(1,P1_i.ndof), zeros(1,DIV_i.ndof),         zeros(1,P1_e.ndof),                   0,                      veci',                   0,                   0; %1
             ];

end