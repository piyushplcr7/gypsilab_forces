addpath(genpath("../../../../"));
clear; 
clc;
format long;
global X;
global W;
load('X3','X');
load('W3','W');

Nvals = 50:200:2000;
sz = size(Nvals,2);
torques_mst = zeros(sz,3); 
torques_bem = zeros(sz,3);
forces_mst = zeros(sz,3);
forces_bem = zeros(sz,3);

for i = 1:sz
    
    % Get the mesh
    [mesh,mesh_in,mesh_out] = tor_tor_mesh(10,3,Nvals(i),10);
    
    % Solve the floating potential problem on mesh
    [Psif,cf] = solve_float_pt_ext(mesh,mesh_in,1,3,'gypsi');
    
    % Definition of FEM spaces and integration rule
    S1_Gamma = fem(mesh,'P1');
    S0_Gamma = fem(mesh,'P0');
    order = 3;
    Gamma = dom(mesh,order);
    
    Op_in = restriction(S0_Gamma,mesh_in);
    
    % Getting the Single Layer matrix V using Gypsilab implementation
    Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); % 0 wave number
    V = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,Gxy,S0_Gamma);
    V = V + 1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'[1/r]',S0_Gamma);
    
    % Getting the Double Layer matrix K using Gypsilab implementation
    GradG = cell(3,1);
    GradG{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
    GradG{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
    GradG{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);
    K = 1/(4*pi)*integral(Gamma,Gamma,S0_Gamma,GradG,ntimes(S0_Gamma));
    K = K +1/(4*pi)*regularize(Gamma,Gamma,S0_Gamma,'grady[1/r]',ntimes(S0_Gamma));
    
    % Defining the mass matrix M
    M = integral(Gamma,S0_Gamma,S0_Gamma);
    
    % Defining the Dirichlet boundary condition
    % Cutoff radius for the BC
    %R = 1; 
    R = 18;
    %assert(R < Rad + s);
    Vin = 1;
    Vout = 0;
    g = @(X) (sqrt(sum(X.^2,2)) > R)* Vout + (sqrt(sum(X.^2,2)) <= R) * Vin;
    
    % Checking the boundary condition
    figure;
    plot(mesh);
    hold on;
    plot(mesh,g(S0_Gamma.dof));
    title('Dirichlet Boundary Condition');
    colorbar;
    
    % Constructing the RHS
    g_N = integral(Gamma,S0_Gamma,g);
    
    % Exterior problem
    Psi = V\((-0.5 * M + K)* (M\g_N));
    
    Psi_in = Op_in * Psi;
    Psi_out = Op_out * Psi;

end
