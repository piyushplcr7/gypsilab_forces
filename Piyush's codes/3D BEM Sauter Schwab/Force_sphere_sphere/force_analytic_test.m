% Testing the force formula

% Testing convergence for some alpha values
alpha1 = acosh(10);
alpha2 = acosh(20);

force1 = [];
force2 = [];

for i = 1 : 500
   force1 = [force1 Gsum(i,alpha1)];
   force2 = [force2 Gsum(i,alpha2)];
end

figure;

plot(1:500,force1);
hold on;
plot(1:500,force2);

%alpha = acosh(1+x);

% Next plot

xvals = 0.001:0.001:0.01;
Gvals = [];
for i = 1 : size(xvals,2)
    Gvals = [Gvals Gsum(100,acosh(1+xvals(i)))];
   
end

figure;

plot(xvals,Gvals);