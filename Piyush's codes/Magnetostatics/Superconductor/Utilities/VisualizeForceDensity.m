% Computing net force using the MST formula and vector potential traces
%
% Inputs
%
% TdA : Dirichlet Trace in the NED Space
% TnA : Neumann Trace in the RWG Space
% Gamma : Integration domain

function out = VisualizeForceDensity(TdA,TnA,Gamma)
    % Constructing the BEM Spaces
    NED = fem(Gamma.msh,'NED');
    RWG = fem(Gamma.msh,'RWG');
    curlNED = NED.curl;
    

    % Quadrature weights and points for integration
    [X,W] = Gamma.qud;
    normals = Gamma.qudNrm;

    % Computing the surface curl of g at the quadrature points
    curlg = reconstruct(TdA,Gamma,curlNED);
    % Computing the Neumann Trace at the quadrature points
    psi = reconstruct(TnA,Gamma,RWG);

    force_n = 0.5 * (curlg .* curlg - dot(psi,psi,2)) .* normals;
    force_t = curlg .* (cross(normals,psi,2));
    
    force = force_n + force_t;
    quiver3wrapper(X,force,'red');

    out = [];
    
end