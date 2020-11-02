%% Computation of int_{0}^2 u^{-alpha}/(u-1) du

alpha = 0.5;
res = 1;

N = 15;

[x,w] = Gauss_Legendre1D(N,0,1);
x1 = x/2 + 3/2;
w1 = w*1/2;
fun1 = @(u)((u.^(-alpha))./(u - 1));

I1 = sum(fun1(x1).*w1);

x2 = pi*x;
w2 = pi*w;
fun2 = @(theta)(1i*(1 + exp(1i*theta)/2).^(-alpha));
I2 = -sum(fun2(x2).*w2);

[x3,w3] = gaussQuadPower(alpha,N);
x3 = x3/2;
w3 = w3.*(1/2).^(1 - alpha);
fun3 = @(u)(1./(u-1));
I3 = sum(fun3(x3).*w3);

I = I1 + I2 + I3 + 1i*pi*res;


matlabGuess = integral(@(u)((u.^(-alpha)-1)./(u-1)),0,2)

%% Computation of int_{0}^b y^{-alpha}/(y-x) dy

N = 10;
alpha = 0.7
X = 0.5
b = 1;

F1 = @(y)((y.^(-alpha) - X.^(-alpha))./(y-X));
CX = X.^(-alpha).*log(abs(b-X)./X);

tic;
matlabGuess = integral(F1,0,b) + CX;
toc;

tic;

switch b < X/2
    case true
        disp('case 1');
        [x,w] = gaussQuadPower(alpha,N);
        x = x.*b./X;
        w = w.*(b./X).^(1-alpha);
        fun = @(u)(1./(u-1));
        I = X.^(-alpha).*sum(fun(x).*w);
    otherwise
        disp('case 2')
        [x,w] = gaussQuadPower(alpha,N);
        x = x.*1/2;
        w = w.*(1/2).^(1-alpha);
        fun1 = @(u)(1./(u-1));
        I1 = X.^(-alpha).*sum(fun1(x).*w);
        
        if b/X < 2
            [x,w] = Gauss_Legendre1D(N,1/2,b/X);
            fun2 = @(u)(((u~=1).*(u.^(-alpha) - 1))./((u~=1).*u-1) + (u==1).*(-alpha));
            I2 = X.^(-alpha).*((sum(fun2(x).*w))+ log(2*abs(b-X)./X));
            I3 = 0;
        else
            disp('case 3');
            [x,w] = Gauss_Legendre1D(N,1/2,2);
            fun2 = @(u)((u.^(-alpha) - 1)./(u-1));
            I2 = X.^(-alpha).*((sum(fun2(x).*w))+ log(2));
            %             if b/X < 5
            %                 % Use a composite rule for the remainder
            %                 [x,w] = compositeGaussLegendre(fix(b/X)+1,N,2,b/X);
            %                 fun3 = @(u)(u.^(-alpha)./(u-1));
            %                 I3 = X.^(-alpha).*sum(fun3(x).*w);
            %             else
            % Use a log scale
            u = log(X);
            fun3 = @(t)(exp(-(alpha)*t)./(1 - exp(u - t)));
            Q = max(fix(2*N*(log(b/(2*X))))+1,N);
            [x,w] = Gauss_Legendre1D(Q,log(2) + u,log(b));
            I3 = sum(fun3(x).*w);
            %             end
        end
        
        
        I = I1 + I2 + I3;
        
end
toc;


disp(I - matlabGuess)
disp(abs((I - matlabGuess)/I));



%% Test the vectorized version:
clear all
close all
clc;
b = 1;
X = 1e-1;
alpha = 0.7;
N = 15;
t1= tic;
I = specialIntegral1(alpha,b,X,N);
t1 = toc(t1);

t2 = tic;
for i = 1:length(X)
    F1 = @(y)((y.^(-alpha) - X(i).^(-alpha))./(y-X(i)));
    CX = X(i).^(-alpha).*log(abs(b-X(i))./X(i));
    
    Imatlab(i,1)  = integral(F1,0,b) + CX;
end
t2 = toc(t2);
disp(max(abs(I - Imatlab)));
plot(X,I);