close all
clear all;
NgaussPhi = 5;
NgaussR = 2;
Ntot = 2*NgaussPhi*NgaussR;
disp(Ntot);
% m = mshDisk(100,1);
% 
% figure; 
% plot(m);
% hold on;
% for i = 1:length(m)
%     T = m.vtx(m.elt(i,:),:);
%     zA = T(1,1) + 1i*T(1,2);
%     zB = T(2,1) + 1i*T(2,2);
%     zC = T(3,1) + 1i*T(3,2);
%     [z,w] = magicQuad(zA,zB,zC,NgaussPhi,NgaussR);
%     plot(z,'g*');
% end
% 

zA = 0.0426 + 1i*(-0.9966);
zB = 0 - 1i*0.8333;
zC = -0.0826 + 1i*(-0.9966);

A = [real(zA);imag(zA)];
A = A/(norm(A)+0.01);
B = [real(zB);imag(zB)];
C = [real(zC);imag(zC)];
C = C/(norm(C)+0.01);

zA = A(1) + 1i*A(2);
zB = B(1) + 1i*B(2);
zC = C(1) + 1i*C(2);

plot([zA zB zC zA],'k');
axis equal
hold on
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k--')
[z,w] = magicQuad(zA,zB,zC,NgaussPhi,NgaussR);
plot(z,'g*');
I = sum(sum(w));

f = @(x,y)(1./sqrt(1 - x.^2 - y.^2));
Icheck = intTri(A,B,C,f);
disp(I - Icheck)
% 
% elt = [1 2 3];
% vtx = [[A; 0]';[B; 0]';[C; 0]'];
% m = msh(vtx,elt);
% Gamma = dom(m,3);
% I= integral(Gamma,@(X)(cos(X(:,1)).*cos(X(:,2))));
% Icheck = intTri(A,B,C,@(x,y)(cos(x).*cos(y)));


