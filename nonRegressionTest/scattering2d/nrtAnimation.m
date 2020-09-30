clean;


[X,Y,Z] = FunR3.XYZ;
theta = rand(1); u = [cos(theta), sin(theta),0]; % direction of propagation.



j = exp(1i*2*pi/3);
p = [1,j,j^2];

m = bnd(mshSquare(100,[1,1]));
x = linspace(-3,3,100);
y = x;
[Xgrid,Ygrid] = meshgrid(x,y);
radiat = [Xgrid(:), Ygrid(:), 0*Ygrid(:)];


Gamma = dom(m,3);
Vh = fem(m,'P1');

k0 = 5;
ks = linspace(k0-3,k0+5,30);
sigma = 1;
amps = exp(-sigma/2*(ks-k0).^2);
figure;
plot(ks,amps);
title('Paquet d''ondes')

tot = cell(length(ks),1);
reg = -1/(2*pi)*regularize(radiat,Gamma,'grady[log(r)]',ntimes(Vh));
for i = 1:length(ks)
    amp = amps(i);
    k = ks(i);
    PW = exp(-1i*k*(u(1)*X + u(2)*Y + u(3)*Z));
    dPW1 = -1i*u(1)*k*PW;
    dPW2 = -1i*u(2)*k*PW;
    dPW3 = -1i*u(3)*k*PW;
    dPW = {dPW1,dPW2,dPW3};
    
    dudn = integral(Gamma,ntimes(Vh),dPW);
    mu = solveHelmholtzNeumann(Gamma,Vh,dudn,k);
    
    dG    = cell(1,3);
    dG{1} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]1',k);
    dG{2} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]2',k);
    dG{3} = @(X,Y) femGreenKernel(X,Y,'grady[H0(kr)]3',k);
    DL = 1i/4*integral(radiat,Gamma,dG,ntimes(Vh));
    
    DL = DL + reg;
    
    tot{i} = amp*(PW(radiat) - DL*mu);
    vals = reshape(tot{i},length(x),length(y));
end


% Animation 
M = radiat;
dt = 0.08;
figure;
hold off;
while true
    vals = 0;
    for i = 1:length(ks)
        tot{i} = tot{i}*exp(-1i*ks(i)*dt);
        vals = vals + tot{i};
    end
    vals = reshape(vals,length(x),length(y));
    imagesc(x,y,real(vals));
    axis xy;
    caxis([-4,4]);
    drawnow;
end




