%% Class Fe
% This class represents spaces of finite element functions. 
% An instance of Fe represents a vector space
% $V_h = $ Span $({\phi_1,...,\phi_n})$ where $\phi_i$ are Finite Element functions
% defined on a mesh. By finite element functions, we mean the following
% conditions:
% 
% - There exists a set of points $x_1,...,x_n$ on the mesh, called dofs,
% such that $\phi_i(x_j) = 0$ for $j \neq i$ and $\phi_i(x_i) = 1$ for all $i$.
% 
% - Each element E of the mesh is the image of the reference domain $S = [0,1]$ by
% an affine mapping $\chi_E$. There exist B (basis) functions
% $\psi_1,...,\psi_B$ defined on S, and points $X_1,...,X_B$ (reference dofs) of S such that
% $\psi_i(X_b) = \delta_{i,b}$. For every $f$ in $V_h$, there holds
% $f = \sum_{E} \sum_{b = 1}^B f(x_{i(E,b)})\phi_b \circ \chi_j^{-1}$ almost
% everywhere.
% 
% - Each dof $x_i$ located on an element E coincides with $\chi_E(X_b)$ for
% some $b$. If $x_i$ is shared by several elements $e_1,...,e_M$, then there
% exists $b_1,...,b_m$ such that $x_i = \chi_{e_m}(X_{b_m})$ for all $m =
% 1..M$.
    
%% Abstract FE and concrete classes
% Fe.m is an abstract class. This means that one cannot create an instance
% of Fe.m. However, one can create concrete classes implementing Fe.m. This
% is due to the fact that many types of finite element functions exist,
% with some specific aspects (number of basis functions, location of the
% dofs,...) but they are involved in Galerkin computations with similar
% algorithms. The algorithms that are common to every type of Fe space are
% written directly in the abstract class. The concrete classes must
% implement only a few function to inherit from Fe. Those are 
% 
% - X,dof_jb = dof(obj), which returns the coordinates of the dofs on the
% mesh and an array dof_jb such that, for each element e_j and index b, 
% dof_jb(j,b) = i where $i$ is such that that $x_i = \chi_{e_j}(X_b)$.
% 
% - r = psi_b(obj,b,x), which, given some reference coordinates x and an
% index b, returns $\psi_b(x)$. 
% 
% - s = name(obj), which simply indicates the name of the finite element
% functions. 


clean;
c = openline(0,2); m = meshCurve(c,2);

% Example 1:
Vh = P0(m); % Space of piecewise constant functions on the mesh. 
[X,dof_jb] = dof(Vh);
disp('P0 dofs :')
disp(X); 
disp('P0 dof_jb')
disp(dof_jb);
figure
plot(m);
hold on
plot(X(:,1),X(:,2),'r*','DisplayName','dofs');
legend show
title('P0 elements')

% Example 2:
Vh = P1(m); % Space of piecewise linear functions on the mesh. 
[X,dof_jb] = dof(Vh);
disp('P1 dofs :')
disp(X); 
disp('P1 dof_jb')
disp(dof_jb);
figure
plot(m);
hold on
plot(X(:,1),X(:,2),'r*','DisplayName','dofs');
legend show
title('P1 elements')

% Example 3:
Vh = P2(m); % Space of piecewise quadratic functions on the mesh. 
[X,dof_jb] = dof(Vh);
disp('P2 dofs :')
disp(X); 
disp('P2 dof_jb')
disp(dof_jb);
figure
plot(m);
hold on
plot(X(:,1),X(:,2),'r*','DisplayName','dofs');
legend show
title('P2 elements');

%% Matrices of some relevant linear maps. 
% In Galerkin Boundary Element methods, the linear systems take the form 
%
% $$ A_{ij} = \int_{D \times D} (K \phi_i)(x) G(x,y) (L\phi_j)(y) \,dx\,dy $$
% 
% where $D$ is a domain of integration and $K$ and $L$ are some linear 
% operators (possibily identity) and $G$ is a possibily singular kernel. Using Gaussian 
% quadratures and ignoring the singularity integrals for the time being, 
% this takes the form
% 
% $$ A_{i,j} \approx \sum_{k,l} \omega_{k,l} (K\phi_i)(x_k) G(x_k, y_l) (L\phi_j)(x_l) $$
% 
% This can be written compactly as $A \approx M_K M_G M_L$ where 
% 
% - $M_G$ is the matrix given by $(G(x_k,y_l))_{k,l}$
% 
% - $M_L$ is the matrix of the linear map which, to the values of a member of
% $u \in V_h$ at the dofs, associates the values of $Lu$ at the points
% $(x_k)$. 
% 
% In the particular case where $L = I_d$, the column $i$ of $M_L$ is simply
% the evalutation of $\phi_i$ at the points $(x_k)$. This matrix can be
% computed by the following command:

% P1 :
Vh = P1(m);
Y = [linspace(0,2,100)', zeros(100,1)];% Query points
X = dof(Vh);
M = dof2Y(Vh,Y);
figure
plot(m);
hold on
plot(X(:,1),X(:,2),'r*')
plot(Y(:,1),full(M));
legend('vertices','dofs','phi_1','phi_2','phi_3')
title('P1 elements')

% P2 :
Vh = P2(m);
Y = [linspace(0,2,600)', zeros(600,1)];% Query points
X = dof(Vh);
M = dof2Y(Vh,Y);
figure
plot(m);
hold on
plot(X(:,1),X(:,2),'r*')
plot(Y(:,1),full(M));
legend('vertices','dofs','phi_1','phi_2','phi_3','phi_4','phi_5')
title('P2 elements')

%%
% If L is not the identity mapping, just compute the matrix by M =
% dof2Y(L(Vh),Y). If $L(\phi_i)$ has several componenents (for example if L
% is defined by $L\phi = \phi \vec{n}$ where $\vec{n}$ is the normal vector
% on the mesh), then M will be returned as a cell with M{k} containing the
% matrix of the k-th component of $L$. 

% Example: tangential gradient of P2 elements
Vh = P2(m);
Y = [linspace(0,2,600)', zeros(600,1)];% Query points
X = dof(Vh);
M = dof2Y(grad(Vh),Y);
figure
plot(X(:,1),X(:,2),'r*')
hold on
plot(Y(:,1),full(M));
legend('dofs','grad(phi_1)','grad(phi_2)','grad(phi_3)','grad(phi_4)','grad(phi_5)')
title('grad(P2) elements')
axis auto
xlim([0,2])
ylim([-4,4])

% Example: multiplication by x of P1 elements
m = meshCurve(circle,7);
Vh = P1(m);
X = dof(Vh);
Y = linspace(m,20);
M = dof2Y(xtimes(Vh),Y);
figure
plot(X(:,1),X(:,2),'r*')
hold on
plot(Y(:,1),full(M{1}));
plot(m)
title('x1 * P1 elements')
figure
plot(X(:,1),X(:,2),'r*')
hold on
plot(Y(:,1),full(M{2}));
plot(m)
title('x2 * P1 elements')






