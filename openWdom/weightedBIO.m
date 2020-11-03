function [Somega,Nomega,rq,loc] = weightedBIO(Gamma,Vh,k,varargin)

p = inputParser;
p.addOptional('lambda',1);
p.addOptional('tol',1e-3);
p.addOptional('Ncompress',15e3);
p.parse(varargin{:});
lambda = p.Results.lambda;
tol = p.Results.tol;
Ncompress = p.Results.Ncompress;

[X,W] = qud(Gamma);
if k == 0
    G = LogKernel(1);
    C = -1/(2*pi);
else
    G = H0Kernel(k);
    C = 1i/4;
end

N = size(X,1);


if N < Ncompress
    Gxy = zeros(N,N);
    for j = 1:N
        Gxy(:,j) = G.eval(norm3D(X-X(j,:)));
    end
else
    a = lambda/N^(2/3);
    [A,rq,loc] = offlineEBD(G,X,X,a,tol);
    Gxy = AbstractMatrix([],@(x)A(x),N,N);
end

% Single layer :
Mu = Vh.uqm(Gamma);
Wx = spdiags(W,0,size(W,1),size(W,1));
S = Mu.'*Wx*Gxy*Wx*Mu;
Sreg = regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);
Somega = C*S  - 1/(2*pi)*Sreg;

% Hypersingular
% Part 1: 
v = omegaDomega(Vh);
Mv = v.uqm(Gamma);
N1 = Mv.'*Wx*Gxy*Wx*Mv;
N1reg = regularize(Gamma,Gamma,v,'[log(r)]',v);
Nomega1 = C*N1 -1/(2*pi)*N1reg;
% Part 2 if k > 0
if k > 0
    w = ntimes(Vh);
    Mw = w.uqm(Gamma);
    omega2 = 1/(Gamma.w)^2;
    omega2X = spdiags(omega2(X),0,N,N);
    N2 = 0;
    for j = 1:3
       N2 = N2 + Mw{j}.'*Wx*omega2X*Gxy*omega2X*Wx*Mw{j}; 
    end
    
    N2reg = regularize(Gamma,Gamma,w,omega2,'omega2[log(r)]',w);
    Nomega2 = C*N2 - 1/(2*pi)*N2reg;
else
    Nomega2 = 0;
end

Nomega = Nomega1 - k^2*Nomega2;


end

