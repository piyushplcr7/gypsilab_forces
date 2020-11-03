clear all;
close all;
clc;


GMRESTOL = 1e-8;
GMRESMAXIT = 500;
TESTMODE = true;
gss = 3;
tolEBD = 1e-3;
    


[X,Y,Z] = FunR3.XYZ;

%% Laplace Segment


c = openline(-1,1);
k = 0;


if TESTMODE
    resID = fopen('resultsSimuNMLaplaceTest.txt','W');
    Ns = [20 80 320 1500];
    lambdas = [0 0 0 10];
    NstartCompress = 1e3;
    
else
    resID = fopen('resultsSimuNMLaplace.txt','W');
    Ns = [200 800 3000 12000];
    lambdas = [0 0 0 1];
    NstartCompress = 5e3;
end

%% 1°) Symm's equation

fprintf(resID,'SYMMS EQUATION \n \n');

for i = 1:length(Ns)
    N = Ns(i);
    fprintf(resID,'N = %s : \n',num2str(N));
    lambdaEBD = lambdas(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    edges = bnd(m);
    % Weight definition :
    omega2 = 1 - X^2;
    omega = sqrt(omega2);
    
    
    
    singVtx = edges.vtx; % Singularities of 1/omega
    singPow = [-1/2;-1/2]; % Power law of the singularities
    sing = {singVtx,singPow};
    
    % Integration domain
    gss = 3;
    Gamma = Wdom(m,gss,1/omega,sing);
    
    
    % Boundary element space
    Vh = P1(m);
    
    k = 0;
    GXY = @(X,Y)femGreenKernel(X,Y,'[log(r)]',k);
    
    % Symm's integral operator
    if N<NstartCompress
        Somega = -1/(2*pi)*(...
            integral(Gamma,Gamma,Vh,GXY,Vh)  ...
            + regularize(Gamma,Gamma,Vh,'[log(r)]',Vh));
    else
        Somega = -1/(2*pi)*(...
            integralEbd(Gamma,Gamma,Vh,'[log(r)]',0,Vh,...
            tolEBD,lambdaEBD)  ...
            + regularize(Gamma,Gamma,Vh,'[log(r)]',Vh));
    end
    
    
    % Mass matrix int_{Gamma} phi_i phi_j/(1-X^2)^{1/2}
    Iomega_1 = integral(Gamma,Vh,Vh);
    
    % Weighted Laplacian
    OmegaDx2 = integral(Gamma,grad(Vh),omega2,grad(Vh));
    
    % Galerkin matrix of the L^2_{1/omega} orthogonal projection on constants
    L = integral(Gamma,0*X + 1);
    C = integral(Gamma,Vh);
    pi0 = @(u)((C'*u)/L*C);
    % Square-root of the Weighted Laplacian
    SQ = @(u)TrefethenSqrtGalerkin(OmegaDx2,5,u,Iomega_1,1,length(Vh)^2);
    
    aux = @(v)(SQ(v) + 2/log(2)*(pi0(v)));
    Prec = @(u)(Iomega_1\(aux(Iomega_1\u)));
    
    % rhs "almost singular"
    rhs = integral(Gamma,Vh,1/sqrt(1/N^2 + X^2));
    %
    t0 = tic;
    [~,~,~,~,RESVEC0]  = gmres(Somega,rhs,[],GMRESTOL,GMRESMAXIT);
    t0 = toc(t0);
    
    t2 = tic;
    [~,~,~,~,RESVEC2]  = gmres(Somega,rhs,[],GMRESTOL,GMRESMAXIT,Prec);
    t2 = toc(t2);
    
    fprintf(resID,'No preconditioner : %s iterations in %s seconds \n',...
        num2str(length(RESVEC0)),num2str(t0));
    fprintf(resID,'Square-root preconditioner : %s iterations in %s seconds \n',...
        num2str(length(RESVEC2)),num2str(t2));
    
end

clear Somega omegaDx2 Iomega_1 rhs;

%% 2°) Hypersingular equation


fprintf(resID,'\n \nHypersingular equation \n\n');

for i = 1:length(Ns)
    
    fprintf(resID,'N = %s : \n',num2str(N));
    N = Ns(i);
    lambdaEBD = lambdas(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    edges = bnd(m);
    % Weight definition :
    [X,Y,Z] = FunR3.XYZ;
    omega2 = 1 - X^2;
    omega = sqrt(omega2);
    
    dOmega{1} = -X./omega;
    dOmega{2} = 0*X;
    dOmega{3} = 0*X;
    
    
    singVtx = edges.vtx; % Singularities of 1/omega
    singPow = [-1/2;-1/2]; % Power law of the singularities
    sing = {singVtx,singPow};
    
    % Integration domain
    Gamma = Wdom(m,gss,1/omega,sing);
    Gamma = Gamma.supplyDw(dOmega);
    
    
    % Boundary element space
    Vh = P1(m);
    
    % Hypersingular operator
    if N < NstartCompress
        Nomega = -1/(2*pi)*(...
            integral(Gamma,Gamma,omegaDomega(Vh),GXY,omegaDomega(Vh))  ...
            + regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh)));
    else
        Nomega = -1/(2*pi)*(...
            integralEbd(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',0,omegaDomega(Vh),...
            tolEBD,lambdaEBD)  ...
            + regularize(Gamma,Gamma,omegaDomega(Vh),'[log(r)]',omegaDomega(Vh)));
        
    end
    
    % Mass matrix int_{Gamma} phi_i phi_j/(1-X^2)^{1/2}
    Iomega = integral(Gamma,Vh,omega2,Vh);
    
    % Weighted Laplacian
    DxOmega2 = integral(Gamma,omegaDomega(Vh),omegaDomega(Vh));
    
    % Square-root of the Weighted Laplacian
    SQ = @(u)TrefethenSqrtGalerkin(DxOmega2,15,u,Iomega,1.5,length(Vh)^2);
    
    
    Prec = @(u)(DxOmega2\(SQ(Iomega\u)));
    
    rhs = integral(Gamma,Vh,omega2*(sqrt((1/N^2 + X^2))));
    
    t0 = tic;
    [X0,~,~,~,RESVEC0]  = gmres(Nomega,rhs,[],GMRESTOL,GMRESMAXIT);
    t0 = toc(t0);
    
    t2 = tic;
    [X2,~,~,~,RESVEC2]  = gmres(Nomega,rhs,[],GMRESTOL,GMRESMAXIT,Prec);
    t2 = toc(t2);
    
    
    fprintf(resID,'No preconditioner : %s iterations in %s seconds \n',...
        num2str(length(RESVEC0)),num2str(t0));
    fprintf(resID,'Square-root preconditioner : %s iterations in %s seconds \n',...
        num2str(length(RESVEC2)),num2str(t2));
    
    
    
end

clear Nomega Iomega DxOmega2 rhs
fclose(resID);



%% Helmholtz segment



if TESTMODE
    resID = fopen('resultsSimuNMHelmholtzSegTest.txt','W');
    ks = [1.25 5 20 80]*pi/2;
    Ns = fix(ks*2*5)+1;
    lambdas = [0 0 10 10];
    NstartCompress = 300;
    
else
    resID = fopen('resultsSimuNMHelmholtzSeg.txt','W');
    ks = [12.5 50 200 800]*pi/2;
    Ns = ks*2*5;
    lambdas = [0 0 0 1];
    NstartCompress = 5e3;
end

for i = 1:length(Ns)
    
    
    
    N = Ns(i);
    k = ks(i);
    fprintf(resID,'\n\nN = %s : \n',num2str(N));
    
    lambdaEBD = lambdas(i);
    m = meshCurve(c,N,'varChange',{@cos,[-pi,0]});
    edges = bnd(m);
    
    % Weight definition :
    
    X1 = edges.vtx(1,:);
    X2 = edges.vtx(2,:);
    % Weight definition :
    [X,Y,Z] = FunR3.XYZ;
    w1 = sqrt((X1(1) - X).^2 + (X1(2) - Y).^2 + (X1(3) - Z).^2);
    w2 = sqrt((X2(1) - X).^2 + (X2(2) - Y).^2 + (X2(3) - Z).^2);
    
    omega2 = w1*w2;
    omega = sqrt(omega2);
    
    dw1{1} = (X - X1(1))./w1;
    dw1{2} = (Y - X1(2))./w1;
    dw1{3} = (Z - X1(3))./w1;
    
    dw2{1} = (X - X2(1))./w2;
    dw2{2} = (Y - X2(2))./w2;
    dw2{3} = (Z - X2(3))./w2;
    
    dOmega = cell(1,3);
    for j = 1:3
        dOmega{j} = (dw1{j}*w2 + dw2{j}*w1)/(2*sqrt(w1*w2));
    end
    
    
    singVtx = edges.vtx; % Singularities of 1/omega
    singPow = [-1/2;-1/2]; % Power law of the singularities
    sing = {singVtx,singPow};
    
    % Integration domain
    Gamma = Wdom(m,gss,1/omega,sing);
    
    Gamma = Gamma.supplyDw(dOmega);
    
    Vh = P1(m);
    
    [Somega,Nomega] = weightedBIO(Gamma,Vh,k,...
        'lambda',lambdaEBD,'tol',tolEBD,'Ncompress',gss*NstartCompress);
    
    
    omegaDx2 = integral(Gamma,grad(Vh),omega2,grad(Vh));
    Omega2 = integral(Gamma,Vh,omega2,Vh);
    Iomega_1 = integral(Gamma,Vh,Vh);
    
    
    K = omegaDx2 - k^2*(Omega2 - Iomega_1);
    keps = k + 1i*0.0025*k^(1/3);
    SQ1 = @(u)(DarbasPadeSqrt(u,20,pi/3,keps,Iomega_1,K));
    PrecSQ1 = @(v)(Iomega_1\SQ1(Iomega_1\v));
    
    
    dxOmega2 = integral(Gamma,omegaDomega(Vh),omegaDomega(Vh));
    Omega2 = integral(Gamma,Vh,omega2^2,Vh);
    Iomega = integral(Gamma,Vh,omega2,Vh);
    
    K = dxOmega2 - k^2*(Omega2 - Iomega);
    D = dxOmega2 - k^2*Omega2;
    SQ2 = @(u)DarbasPadeSqrt(u,20,pi/3,keps,Iomega,K);
    PrecSQ2 = @(u)(D\SQ2(Iomega\u));
    
    PrecSN = @(u)(Iomega_1\(Somega*(Iomega\u)));
    PrecNS = @(u)(Iomega\(Nomega*(Iomega_1\u)));
    
    theta_inc = 0;
    PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
    rhs1 = -integral(Gamma,Vh,PW);
    
    
    fprintf(resID,'\n\nDirichlet problem \n\n');
    
    t0 = tic;
    [~,~,~,~,RESVEC0] = gmres(Somega,rhs1,[],GMRESTOL,GMRESMAXIT);
    t0 = toc(t0);
    t2 = tic;
    [~,~,~,~,RESVEC2] = gmres(Somega,rhs1,[],GMRESTOL,GMRESMAXIT,PrecSQ1);
    t2 = toc(t2);
    t3 = tic;
    [~,~,~,~,RESVEC3] = gmres(Somega,rhs1,[],GMRESTOL,GMRESMAXIT,PrecNS);
    t3 = toc(t3);
    
    fprintf(resID,'No prec. : %s iterations in %s seconds \n',...
        num2str(length(RESVEC0)),num2str(t0));
    fprintf(resID,'Square-root : %s iterations in %s seconds \n',...
        num2str(length(RESVEC2)),num2str(t2));
    fprintf(resID,'Calderon : %s iterations in %s seconds \n',...
        num2str(length(RESVEC3)),num2str(t3));
   

    fprintf(resID,'\n\nNeumann problem \n\n');
   
    
    theta_inc = pi/4;
    PW = exp(1i*k*(X*cos(theta_inc) + Y*sin(theta_inc)));
    omega2dPW{1} = omega2*1i*k*cos(theta_inc)*PW;
    omega2dPW{2} = omega2*1i*k*sin(theta_inc)*PW;
    omega2dPW{3} = 0*PW;
    
    rhs2 = integral(Gamma,ntimes(Vh),omega2dPW);
    
    t0 = tic;
    [~,~,~,~,RESVEC0] = gmres(Nomega,rhs2,[],GMRESTOL,GMRESMAXIT);
    t0 = toc(t0);
    t2 = tic;
    [~,~,~,~,RESVEC2] = gmres(Nomega,rhs2,[],GMRESTOL,GMRESMAXIT,PrecSQ2);
    t2 = toc(t2);
    t3 = tic;
    [~,~,~,~,RESVEC3] = gmres(Nomega,rhs2,[],GMRESTOL,GMRESMAXIT,PrecSN);
    t3 = toc(t3);
    
    
    fprintf(resID,'No prec. : %s iterations in %s seconds \n',...
        num2str(length(RESVEC0)),num2str(t0));
    fprintf(resID,'Square-root : %s iterations in %s seconds \n',...
        num2str(length(RESVEC2)),num2str(t2));
    fprintf(resID,'Calderon : %s iterations in %s seconds \n',...
        num2str(length(RESVEC3)),num2str(t3));
    
end


fclose(resID);






