function Y = forces(mesh)

figure;
title('Volume mesh');
plot(mesh);
figure;
bnd_mesh = mesh.bnd;
normals = bnd_mesh.nrm;
centers = bnd_mesh.ctr;
title('Boundary mesh');
plot(bnd_mesh);
hold on;

%% Creating projection matrix linking boundary triangles to boundary edges

% Number of elements in boundary and volume mesh
Nelt_b = size(bnd_mesh.elt,1);
Nelt_v = size(mesh.elt,1);
tr_opr = sparse(Nelt_b, Nelt_v);

% Looping over all boundary elements and finding the parent volume element
for i = 1 : Nelt_b
    bvindx = bnd_mesh.elt(i,:);
    % Nodes of the boundary element
    bvpts = bnd_mesh.vtx(bvindx,:);
    % Checking which volume mesh nodes match with selected boundary nodes
    distpt1 = vecnorm(mesh.vtx-ones(size(mesh.vtx,1),1).*bvpts(1,:),2,2);
    distpt2 = vecnorm(mesh.vtx-ones(size(mesh.vtx,1),1).*bvpts(2,:),2,2);
    % Indices for matching nodes
    j1 = find(distpt1 == 0);
    j2 = find(distpt2 == 0);
    % Checking which volume element contains both nodes
    for j = 1:Nelt_v
       cur_elt = mesh.elt(j,:);
       % Checking if coincides
       check1 = cur_elt(1) == j1 | cur_elt(2) == j1 | cur_elt(3) == j1;
       check2 = cur_elt(1) == j2 | cur_elt(2) ==j2 | cur_elt(3) == j2;
       if check1 & check2
           disp('Found match');
          tr_opr(i,j)=1; 
       end
    end
end
quiver(centers(:,1),centers(:,2),normals(:,1),normals(:,2));

%% Setting up the discrete system

% Quadrature order
order = 3;

% Creating the integration domain
Omega = dom(mesh, order);
Gamma = dom(bnd_mesh, order);

% FEM space without the dirichlet boundary condition
S1_Omega = fem(mesh,'P1');
% FEM space with the homogeneous Dirichlet BC
S10_Omega = dirichlet(S1_Omega, bnd_mesh);
S1_Gamma = fem(bnd_mesh,'P1');

% LHS matrix
A = integral(Omega,grad(S10_Omega),grad(S10_Omega));

% Restriction operator for volume dofs to boundary dofs
P = restriction(S1_Omega,bnd_mesh); % Gives a matrix with size S1_Gamma.Ndof X S1_Omega.Ndof

% Getting the RHS vector
g_N = g(S1_Gamma.dof);
G_N = P' * g_N; % Extension by zero

v = -integral(Omega,grad(S10_Omega),grad(S1_Omega))* G_N;

% This is the offset solution u* 
sol = A\v; 
% Sol contains N_space - N_bnd elements, how to get a solution of size N_space?
% Using Martin's approach of "elimination"
Q = elimination(S1_Omega,bnd_mesh); % Matrix of size S1_Omega.Ndof X (S1_Omega.Ndof-S1_Gamma.Ndof)

% Solution with full size N_space, extension by zero
sol_ex = Q * sol;

% Obtaining the solution u by undoing the offset
u = sol_ex + G_N;

%% Plotting the solution
pts = S1_Omega.dof;
X = pts(:,1);
Y = pts(:,2);
figure;
scatter3(X,Y,u);
figure;
r = vecnorm(pts,2,2);
u_exact = log(r/1.5)/log(0.5/1.5);
scatter3(X,Y,abs(u_exact-u));
title('Solution error');
figure;
graph(S1_Omega,u);
title('Electrostatic potential');

%% Computing gradu using the volume solution u

% Using Martin's trick to solve for gradu using an integral equation
S0_Omega = fem(mesh,'P0');
Omega = dom(mesh,3);
gradu = cell(3,1);
% Computing the LHS Mass matrix
M00 = integral(Omega,S0_Omega,S0_Omega);
% Computing the RHS Matrix X vector

% Getting coefficients of gradu in terms of basis functions of gradspace
gradu{1} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,1)) * u);
gradu{2} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,2)) * u);
gradu{3} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,3)) * u);

% Visualizing the electric field
figure;
plot(mesh);
title('Electric field');
hold on;
%centers = mesh.ctr;
centers = S0_Omega.dof;
quiver(centers(:,1),centers(:,2),gradu{1},gradu{2});

%% Computing gradu.n, need to restrict P0 from volume to P0 on the boundary

% Natural fem space for gradu on the boundary: P0
S0_Gamma = fem(bnd_mesh,'P0'); % S0_Gamma
pts = S0_Gamma.dof;
% Using Martin's function traceP0P1
% Trace operator from S0_Omega to S0_Gamma`
%[tr_op,~] = trace_P0P1(S0_Omega);

% 
tr_op = tr_opr;

bgradu = cell(3,1);

bgradu{1} = tr_op * gradu{1};
bgradu{2} = tr_op * gradu{2};
bgradu{3} = tr_op * gradu{3};
quiver(pts(:,1),pts(:,2),bgradu{1},bgradu{2},'red');
bnd_norms = bnd_mesh.nrm;

% Computing the Neumann trace
Neu_tr = bgradu{1}.*bnd_norms(:,1) + bgradu{2}.*bnd_norms(:,2);

% Plotting the Neumann trace
figure;
plot(bnd_mesh);
title('Neumann Trace');
hold on;
pts = S0_Gamma.dof;
scatter3(pts(:,1),pts(:,2),Neu_tr);

%% Computing the volume based expression for force
% Getting quadrature weights and nodes
[QP,QW,elmat] = Omega.qud;
% Getting the unknown to quadrature matrix
uqmat = S0_Omega.uqm(Omega);

% Need to define the function Psi \in P1 and then get coefficients of 
% grad Psi \in P0

% Extracting the elements corresponding to the inner boundary
% element centers are obviously ordered according to the elements
bnd_centers = mesh.bnd.ctr;
dist_bnd_centers = vecnorm(bnd_centers,2,2);
inner_elt_idx = find( dist_bnd_centers < 1.2 );
% mesh.sub(idx) creates a mesh containing elements specified by
% mesh.elt(idx)
inner_bnd_mesh = bnd_mesh.sub(inner_elt_idx);
S1_Gamma_in = fem(inner_bnd_mesh,'P1');

% Correct way to find out the inner indices
%bnd_centers = mesh.bnd.ctr;
%inner_idx = find(abs(bnd_centers(:,1))<1.99 & (abs(bnd_centers(:,2)) < 1.99));

% Creating a restriction operator, size S1_Gamma_in.ndof X S1_Omega.ndof
% Dirty trick for cut off function
%P = restriction(S1_Omega,inner_bnd_mesh);
%Psi = P' * ones(S1_Gamma_in.ndof,1);

% Smooth cut off function using cos^2
S1_Omega_dofs = S1_Omega.dof;
Psi = cutoff(S1_Omega_dofs(:,1),S1_Omega_dofs(:,2));

% Plotting this function Psi
%figure;
%graph(S1_Omega,Psi)

% Using the same method used for finding the gradient of u:
Omega = dom(mesh,3); % Trying quadrature order 3 this time, might be useless
M00 = integral(Omega,S0_Omega,S0_Omega);
Psi_grad = cell(3,1);
Psi_grad{1} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,1)) * Psi);
Psi_grad{2} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,2)) * Psi);
Psi_grad{3} = M00\(integral(Omega,S0_Omega,grad(S1_Omega,3)) * Psi);

Mat_Psi_grad = [Psi_grad{1}, Psi_grad{2}, Psi_grad{3}];
Mat_u_grad = [gradu{1}, gradu{2}, gradu{3}];

% Matrices storing the evaluations of these gradient functions at
% quadrature points
Psi_grad_at_QP = uqmat * Mat_Psi_grad;
u_grad_at_QP = uqmat * Mat_u_grad;

% Finally evaluating the volume based force formula
f_vol = cell(3,1);
f_vol{1} = sum(QW .* sum(Psi_grad_at_QP .* u_grad_at_QP,2) .* u_grad_at_QP(:,1)) - 0.5 * sum(QW .* sum(u_grad_at_QP .* u_grad_at_QP,2) .* Psi_grad_at_QP(:,1));
f_vol{2} = sum(QW .* sum(Psi_grad_at_QP .* u_grad_at_QP,2) .* u_grad_at_QP(:,2)) - 0.5 * sum(QW .* sum(u_grad_at_QP .* u_grad_at_QP,2) .* Psi_grad_at_QP(:,2));
f_vol{3} = sum(QW .* sum(Psi_grad_at_QP .* u_grad_at_QP,2) .* u_grad_at_QP(:,3)) - 0.5 * sum(QW .* sum(u_grad_at_QP .* u_grad_at_QP,2) .* Psi_grad_at_QP(:,3));

%% Computing the integral 1/2 \int_{\Gamma}(gradu.n)^2 n dS

% Neumann trace lies in S0_Gamma, Neu_tr are the corresponding
% coefficients
Gamma_in = dom(inner_bnd_mesh,1);
S0_Gamma_in = fem(inner_bnd_mesh,'P0');
Mbin_00 = integral(Gamma_in,S0_Gamma_in,ntimes(S0_Gamma_in));
% Restricting the Neumann trace to S0_Gamma_in
S0_Gamma = fem(bnd_mesh,'P0');
% Size S0_Gamma_in_Ndof X S0_Gamma_Ndof
opr = restriction(S0_Gamma,inner_bnd_mesh);
Neu_in = opr * Neu_tr;
F_bnd = zeros(3,1);
F_bnd(1) = 0.5 * Neu_in' * Mbin_00{1} * Neu_in;
F_bnd(2) = 0.5 * Neu_in' * Mbin_00{2} * Neu_in;
F_bnd(3) = 0.5 * Neu_in' * Mbin_00{3} * Neu_in;

%% Plotting the new cut off function along with the mesh
figure;
plot(mesh);
hold on;
x = -3:0.05:3;
y=x;
[X,Y] = meshgrid(x,y);
surf(X,Y,cutoff(X,Y));

%% Saving results to a file
fileID = fopen('results.txt','a');
fprintf(fileID,'%d, %d; %d, %d \n',F_bnd(1),F_bnd(2),f_vol{1},f_vol{2});
fclose(fileID);
end

%% Function describing the boundary conditions

function Y = g(X)
    g_o = 0; % Outer boundary condition
    g_i = 1; % Inner boundary condition
    logical_outer = vecnorm(X,2,2)>1.2;
    Y = g_o * logical_outer + g_i * not(logical_outer);
    %Y = X(:,1) + X(:,2);
end

function y = cutoff(x,y)
    r1 = 1.2;
    r2 = 1.9; 
    r = sqrt(x.^2+y.^2);
    y = (r<=r1)*1 + (r>=r2)*0 + (r>r1 & r<r2).*cos((r-r1)/(r2-r1)*pi/2).^2;
end















