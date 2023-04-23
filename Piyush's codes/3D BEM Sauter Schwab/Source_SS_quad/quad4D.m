function [Xss, Wss] = quad4D(nq)

if nq == 2

x = [0.5*(1-1/sqrt(3)); 0.5*(1+1/sqrt(3))];
w = [0.5 ; 0.5];

elseif nq == 3
% N = 3 % ACCURACY 1E-16 CONSTANT KERNEL
x = [0.5*(1-sqrt(3/5)) ; 0.5 ; 0.5*(1+sqrt(3/5))];
w = [5/18 ; 4/9 ; 5/18];
% % 

elseif nq == 4
a = sqrt(3/7 + 2/7*sqrt(6/5));
b = sqrt(3/7 - 2/7*sqrt(6/5));
w1 = (18-sqrt(30))/72;
w2 = (18+sqrt(30))/72;
x = [0.5*(1-a) ; 0.5*(1-b) ; 0.5*(1+b) ; 0.5*(1+a)] ;
w = [w1 ; w2 ; w2 ; w1];

elseif nq == 5
a = 1/3*sqrt(5 + 2*sqrt(10/7));
b = 1/3*sqrt(5 - 2*sqrt(10/7));
w1 = (322-13*sqrt(70))/1800;
w2 = (322+13*sqrt(70))/1800;
x = [0.5*(1-a) ; 0.5*(1-b) ; 0.5 ; 0.5*(1+b) ; 0.5*(1+a)] ;
w = [w1 ; w2 ;64/225 ; w2 ; w1];

elseif nq == 6
%% N = 6
wy = [0.3607615730481386 	0.6612093864662645;
0.3607615730481386 	-0.6612093864662645;
0.4679139345726910 	-0.2386191860831969;
0.4679139345726910 	0.2386191860831969;
0.1713244923791704 	-0.9324695142031521;
0.1713244923791704 	0.9324695142031521];

wgss = wy(:, 1); wgss = wgss / 2.0; w = wgss;
x = wy(:, 2); 
x = 0.5 * (1 + x); 
end


% nq = size(w, 1);

X = zeros(nq^4, 4);
W = zeros(nq^4, 1);

count = 1;
for i1 = 1:nq
   for i2 = 1:nq
      for i3 = 1:nq
         for i4 = 1:nq
          X(count, :) = x([i1 i2 i3 i4])';
          W(count)    = prod(w([i1 i2 i3 i4]));
          count = count + 1;
         end
      end
    end
end
   

[Xss{1}, Wss{1}] = quad1_bndry(X, W); % identical

[Xss{2}, Wss{2}] = quad2_bndry(X, W); % common edge

[Xss{3}, Wss{3}] = quad3_bndry(X, W); % common vertex

[Xss{4}, Wss{4}] = quad4_bndry(X, W); % far away


end