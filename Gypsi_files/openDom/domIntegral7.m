function I = domIntegral7(data)
%-------------------------------------------------------------------------%
% int_{Omega x Omega} v(x)G(x,y)w(y)dydx
% with H-matrix compression
% and correction matrix for near field interactions.
% I = V(G - C)W*


% Domain with quadrature
Xdom   = data{1};
[X,Wx] = Xdom.qud;
Nx     = size(X,1);
Wx     = spdiags(Wx,0,Nx,Nx);

% Domain with quadrature
Ydom   = data{2};
[Y,Wy] = Ydom.qud;
Ny     = size(Y,1);
Wy     = spdiags(Wy,0,Ny,Ny);

% Finite element matrix with integration
u  = data{3};
Mu = u.uqm(Xdom);
if iscell(Mu)
    Mu{1} = Mu{1}' * Wx;
    Mu{2} = Mu{2}' * Wx;
    Mu{3} = Mu{3}' * Wx;
else
    Mu = Mu' * Wx;
end

% Green kernel
green = data{4};

% Finite element matrix with integration
v  = data{5};
Mv = v.uqm(Ydom);
if iscell(Mv)
    Mv{1} = Wy * Mv{1};
    Mv{2} = Wy * Mv{2};
    Mv{3} = Wy * Mv{3};
else
    Mv = Wy * Mv;
end

% Accuracy
tol = data{6};

% Correction matrix: Singular integrals
corr = data{7};

% H-Matrix Integration
if iscell(Mu) && ~iscell(green) && ~iscell(Mv)
    I{1} = hmx(u.unk,v.unk,Mu{1},X,green,Y,Mv,tol,corr);
    I{2} = hmx(u.unk,v.unk,Mu{2},X,green,Y,Mv,tol,corr);
    I{3} = hmx(u.unk,v.unk,Mu{3},X,green,Y,Mv,tol,corr);

elseif ~iscell(Mu) && iscell(green) && ~iscell(Mv)
    I{1} = hmx(u.unk,v.unk,Mu,X,green{1},Y,Mv,tol,corr);
    I{2} = hmx(u.unk,v.unk,Mu,X,green{2},Y,Mv,tol,corr);
    I{3} = hmx(u.unk,v.unk,Mu,X,green{3},Y,Mv,tol,corr);

elseif ~iscell(Mu) && ~iscell(green) && iscell(Mv)
    I{1} = hmx(u.unk,v.unk,Mu,X,green,Y,Mv{1},tol,corr);
    I{2} = hmx(u.unk,v.unk,Mu,X,green,Y,Mv{2},tol,corr);
    I{3} = hmx(u.unk,v.unk,Mu,X,green,Y,Mv{3},tol,corr);

else
    I = hmx(u.unk,v.unk,Mu,X,green,Y,Mv,tol,corr);
end


end