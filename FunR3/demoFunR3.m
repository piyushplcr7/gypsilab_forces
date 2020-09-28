%% class FunR3
% 
% Class representing the algebra of functions 
%
% $$f : R^2 \to C$$
%
% and implementing its operation by overriding the operators +, -, *,... of
% the native language. The purpose of this class is also to provide a type 
% for this kind of object. 

f1 = FunR3; % Identical 0 function. 
f2 = FunR3(1); % constant function equal to 1. 
f3 = FunR3(@(Z)(cos(Z(:,1) + 5*sin(10*Z(:,2))))); % Construction by function handle
% One can provide an expression for f3:
% This expression is then used to assign expressions to other instances
% constructed from f3.
f4 = 2*f3 + f2; 
disp(f4);

X = FunR3.X; Y = FunR3.Y; % functions (X,Y) -> X, (X,Y) -> Y.
% FunR3 also supports conjugation and composition with elementary
% functions:
PW = conj(exp(1i*(5*X + 4*Y))/(1 + log(2*Y^2)^2)); 