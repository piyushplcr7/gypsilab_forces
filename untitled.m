vtx = [0 0 0;0 1 0; 1 1 0; 1 0 0];
elt = [1 2 3; 1 3 4];

% vtx = [0 0 0;1 0 0; 0 1 0];
% elt = [1 2 3];


m = msh(vtx,elt);

Gamma = dom(m,7);
Vh = P1(m);

Mh = integral(Gamma,Vh,Vh);
Kh = integral(Gamma,grad(Vh),grad(Vh)) + integral(Gamma,Vh,Vh);

[P,D] = eig(full(Mh\Kh));

