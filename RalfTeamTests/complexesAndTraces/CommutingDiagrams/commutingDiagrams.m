%% Commuting diagrams:

% The following diagram is commuting:


% S^{1,0}(Omega) -> Ned(Omega)   ->   RT(Omega) -> S^{0,-1}(Omega)
%               grad           curl            div
%     |                 |               |               
%     |gamma     1      | piT     2     | pi_n        
%     |                 |               |              
%     V                 V               V               
% S^{1,0}(Gamma) -> Ned(Gamma) -> S^{0,-1}(Gamma) 
%     |         grad    |      div(nx)  |           
%     |Id        3      | nx     4      | Id
%     V                 V               V
% S^{1,0}(Gamma) -> RT(Gamma) -> S^{0,-1}(Gamma)
%               nxgrad       div               

% In the sense that A o B = C o D whenever A and C are two arrows starting
% from the same space and B and D are two arrows pointing to the same
% space.

% When the topology of Omega is trivial, the horizontal seuquences are
% exact. A sequence
% A -> B -> C
%   f    g
% is exact if, for b in B, g(b) = 0 <=> b = f(a) for some a in A. 
% Diagrams 3 and 4 are obvious!


%% Diagram 1:


m = mshCube(30,[1 1 1]);
dm = m.bnd;
Vh = fem(m,'P1');
Wh = fem(dm,'P1');
Xh = fem(m,'NED');

gamma = trace_P0P1(Vh);
piT = tangentialTrace_NED(Xh);
grad3D = grad_P1(Vh);
grad2D = grad_P1(Wh);

Comm = piT*grad3D - grad2D*gamma;
disp(max(max(abs(Comm))))

%% Diagram 2:

Vh = fem(m,'NED');
Wh = fem(dm,'NED');
Xh = fem(m,'RWG');

piT = tangentialTrace_NED(Vh);
pin = normalTrace_RT(Xh);
curl3D = curl_NED(Vh);
curl2D = divnx_NED(Wh);

Comm = pin*curl3D - curl2D*piT;
disp(max(max(abs(Comm))))

%% Diagram 3:

Vh = fem(dm,'P1');
Wh = fem(dm,'P1');
Xh = fem(dm,'NED');

I = fem_eye(Vh);
ncross = nx_RWGNED(Xh);
grad2D = grad_P1(Vh);
nxgrad2D = nxgrad_P1(Wh);

Comm = ncross*grad2D - nxgrad2D*I;
disp(max(max(abs(Comm))))

%% Diagram 4:

Vh = fem(dm,'NED');
Wh = fem(dm,'RWG');
Xh = fem(dm,'P0');

ncross = nx_RWGNED(Vh);
I = fem_eye(Xh);
divncross = divnx_NED(Vh);
div2D = div_RT(Wh);

Comm = I*divncross - div2D*ncross;
disp(max(max(abs(Comm))))

