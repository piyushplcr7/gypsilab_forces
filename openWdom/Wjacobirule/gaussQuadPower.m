function [x,w] = gaussQuadPower(alpha,n,B)

if nargin == 2
    A = 0; B = 1;
else
   A = 0;
   if abs(B) < 1e-13
       x = zeros(n,1);
       w = 0*x;
       return;
   end
   if B < 0
       A = B;
       B = 0;
   end
   
end

[x,w] = jacobi_rule(n,0,-alpha,A,B);
if A<0
    x = A - flipud(x); 
    w = flipud(w);
end


end
