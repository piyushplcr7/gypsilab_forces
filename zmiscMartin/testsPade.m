% Tests Padé. 

close all;
clear all;
clc;
theta = pi/2;
Np = 15;


[ C0,Aj,Bj ] = rotatingPadeRacine(Np,theta);

m = mshDisk(5000,5);
X = m.vtx(:,1)+1i*m.vtx(:,2);
rX = 0*X;
for i = 1:length(X)
    rX(i) = C0;
    for j = 1:length(Aj)
        rX(i) = rX(i) + Aj(j)*(X(i)-1)/(1 + Bj(j)*(X(i)-1));
    end
end

err = abs(rX - sqrt(X));
m.vtx(:,3) = 20*log(err + eps);
plotOn(m,20*log(err + eps));

X = linspace(-1000,1000,200);
X1 = (X+0.0000i)/50;
rX1 = 0*X;
for i = 1:length(X)
    rX1(i) = C0;
    for j = 1:length(Aj)
        rX1(i) = rX1(i) + Aj(j)*(X1(i)-1)/(1 + Bj(j)*(X1(i)-1));
    end
end
err1 = abs(sqrt(50)*rX1 - sqrt(X));
figure;
plot(X,err1);

X = linspace(-100,100,10000);
X1 = -X-1+0.0101i;
rX1 = 0*X1;

for i = 1:length(X)
    rX1(i) = C0;
    for j = 1:length(Aj)
        rX1(i) = rX1(i) + Aj(j)*X1(i)/(1 + Bj(j)*X1(i));
    end
end
err1 = abs(1i*rX1 - sqrt(X));
figure;
plot(X,err1);