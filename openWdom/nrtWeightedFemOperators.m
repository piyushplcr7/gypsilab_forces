

m = mshCircle(500,1);
c = m.ctr;
m = m.sub(c(:,2) > 0 );
proj = @(X)([X(:,1) 0*X(:,2) 0*X(:,3)]);
m = fct(m,proj);
figure; 
plot(m); 
m = swap(m);
title('Adapted mesh with cosine change of variable');



[X,Y,Z] = FunR3.XYZ; % creates the handles [X,Y,Z] -> X, [X,Y,Z] -> Y and [X,Y,Z] -> Z
omega = sqrt(1 - X^2); % FunR3 is able to interpret this expression
singVtx = [[-1,0,0];[1,0,0]]; % Singularities of 1/omega
singPow = [-1/2;-1/2]; % Power law of the singularities
sing = {singVtx,singPow};
gss = 7; 
Gamma = Wdom(m,gss,1/omega,sing); % creates a wDom < dom object. 
Vh = P2(m); % Creates a Fe < fem object. This should change later.
Iomega_1 = integral(Gamma,Vh,Vh);
omegaDx2 = integral(Gamma,grad(Vh),omega^2,grad(Vh));

% the generalized eigenvalue problem int omega dxu dxv = lambda int uv/omega
% has solutions given by lambda_n = n^2, and u_n = Tn 
[~,D] = eig(full(Iomega_1\omegaDx2));
d = sort(diag(D),'ascend');
disp(d(1:10))
err = norm(d(1:10)-((0:9).^2)',2);
fprintf('\n Err = %s \n\n',num2str(err));

Iomega = integral(Gamma,Vh,omega^2,Vh);

dxOmega2 = integral(Gamma,grad(Vh),omega^4,grad(Vh)) ...
    + integral(Gamma,xtimes(Vh),xtimes(Vh)) ...
    - integral(Gamma,xtimes(Vh),omega^2,grad(Vh))...
    - integral(Gamma,grad(Vh),omega^2,xtimes(Vh));

% the generalized eigenvalue problem int 1/omega * (omega dx omega) u (omega dx omega) v 
% = lambda int omega uv 
% has solutions given by lambda_n = (n+1)^2, and u_n = Un 
[P,D] = eig(full(Iomega\dxOmega2));
d = sort(diag(D),'ascend');
disp(d(1:10))
err = norm(d(1:10)-((1:10).^2)',2);
fprintf('\n Err = %s \n\n',num2str(err));


