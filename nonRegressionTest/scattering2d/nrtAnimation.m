clean;


[X,Y,Z] = FunR3.XYZ; 
theta = rand(1); u = [cos(theta), sin(theta),0]; % direction of propagation.

j = exp(1i*2*pi/3);
p = [1 j j^2];
m = polymesh(p,20);
plot(m);
axis equal;
title('Obstacle')

Gamma = dom(m,3);
Vh = fem(m,'P1');

k0 = 10;
ks = linspace(k0 - 1.5,k0+1.5,20);
sigma = 5;
amps = exp(-sigma/2*(ks-k0).^2);
figure;
plot(ks,amps);
title('Paquet d''ondes')

for i = 1:length(ks)
    amp = amps(i);
    k = ks(i);
    PW = exp(1i*(u(1)*X + u(2)*Y + u(3)*Z));
    dPW1 = 1i*u(1)*PW;
    dPW2 = 1i*u(2)*PW;
    dPW3 = 1i*u(3)*PW;
    dPW = {dPW1,dPW2,dPW3};
    
    L1 = integral(Gamma,grad(Vh),dPW);
    L2 = integral(Gamma,grad(P1(m)),dPW);
    
end