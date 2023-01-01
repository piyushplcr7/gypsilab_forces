% Computing net force on boundary Gamma using the MST formula
% 
% Gamma: Dom object representing the boundary of the integration domain
% TD: Coefficients representing the dirichlet trace of the vector potential
% TN: Coefficients representing the Neumann trace of the vector potential
%
function F = magnetostatics_forces(mesh,TD,TN,mu)
    % Dom object for integration
    Gamma = dom(mesh,3);

    % Relevant BEM spaces
    NED = fem(mesh,'NED'); % Dirichlet trace space
    RWG = fem(mesh,'RWG'); % Neumann trace space

    % Integrating the first term, g is the Dirichlet trace
    % 1/(2 mu) \int_{Gamma} (curl_Gamma \vec{g})^2 \vec{n} dS 
    scurlNED = NED.curl;
    uqmat = scurlNED.uqm(Gamma);
    [X,WX,~] = Gamma.qud;
    NX = size(X,1);
    WXmat = spdiags(WX,0,NX,NX); 
    mat1 = uqmat' * WXmat;
    normals = Gamma.qudNrm;
    mat2_1 = normals(:,1).*uqmat;
    mat2_2 = normals(:,2).*uqmat;
    mat2_3 = normals(:,3).*uqmat;

    Mf1 = mat1*mat2_1;
    Mf2 = mat1*mat2_2;
    Mf3 = mat1*mat2_3;

    %MNN = integral(Gamma,NED.curl,NED.curl);
    F1 = 1/2/mu * [ TD' * Mf1 * TD;...
                    TD' * Mf2 * TD;...
                    TD' * Mf3 * TD;];

    % Integrating the second term, lambda is the Neumann trace
    % -mu/2 \int_{Gamma} || \vec{\lambda}||^2 \vec{n} dS
    uqmat_RWG = RWG.uqm(Gamma);
    F2 = [0;0;0];
    for i = 1:3
        uqmat = uqmat_RWG{i};
        mat1 = uqmat' * WXmat;
        mat2_1 = normals(:,1).*uqmat;
        mat2_2 = normals(:,2).*uqmat;
        mat2_3 = normals(:,3).*uqmat;

        Mf1 = mat1*mat2_1;
        Mf2 = mat1*mat2_2;
        Mf3 = mat1*mat2_3;

        F2 = F2 - mu/2 * [TN' * Mf1 * TN;...
                          TN' * Mf2 * TN;...
                          TN' * Mf3 * TN]; 
    end

    % Integrating the third term
    % \int_{Gamma} curl_{Gamma}\vec{g} (\vec{n} x \vec{lambda}) dS
    Mcell = integral(Gamma,scurlNED,NED);
    F3 = -[TD' * Mcell{1} * TN;...
          TD' * Mcell{2} * TN;...
          TD' * Mcell{3} * TN];


    F = F1+F2+F3;
end