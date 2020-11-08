function [Somega,Nomega,rq,loc] = weightedBIO(Gamma,Vh,k,varargin)

p = inputParser;
p.addOptional('lambda',1);
p.addOptional('tol',1e-3);
p.addOptional('Ncompress',100);
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
    rq = 0;
    loc = 0;
else
    a = lambda/N^(2/3);
    disp('Computing EBD')
    [mv,rq,loc,~] = offlineEBD(G,X,X,a,tol);
    Gxy = AbstractMatrix([],mv,N,N);
end

% Single layer :
Mu = Vh.uqm(Gamma);
Wx = spdiags(W,0,size(W,1),size(W,1));
MWu = Wx*Mu;
S = MWu'*(Gxy*MWu);
disp('Regularization SL')
Sreg = regularize(Gamma,Gamma,Vh,'[log(r)]',Vh);
Somega = C*S  - 1/(2*pi)*Sreg;

% Hypersingular
% Part 1: 
v = omegaDomega(Vh);
Mv = v.uqm(Gamma);
MWv = Wx*Mv;
N1 = MWv'*(Gxy*MWv);
disp('Regularization N1')
N1reg = regularize(Gamma,Gamma,v,'[log(r)]',v);
Nomega1 = C*N1 -1/(2*pi)*N1reg;
% Part 2 if k > 0
if k > 0
    disp('Regularization N2')
    w = ntimes(Vh);
    Mw = w.uqm(Gamma);
    omega2 = 1/(Gamma.w)^2;
    Womega2X = spdiags(W.*omega2(X),0,N,N);
    N2 = 0;
    for j = 1:2
       MWj = Womega2X*Mw{j};
       N2 = N2 + MWj.'*(Gxy*MWj); 
    end
    N2reg = regularize(Gamma,Gamma,w,omega2,'omega2[log(r)]',w);
    Nomega2 = C*N2 - 1/(2*pi)*N2reg;
else
    Nomega2 = 0;
end

Nomega = Nomega1 - k^2*Nomega2;


end

