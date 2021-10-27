clear all
close all
theta = 0:0.01:2*pi;
plot(0.3+0.35*cos(theta)+0.1625*cos(2*theta),0.5+0.35*sin(theta),'k');
hold on
axis('equal')
plot([-2,2],[2,2],'k')
plot([-2,2],[-2,-2],'k')
xlim([-2.5,2.5])
ylim([-2.5,2.5])
plot([2,2],[-2,2],'k')
plot([-2,-2],[-2,2],'k')
t=text(-1,-1,'$\Omega$','Interpreter','latex')
t.FontSize = 14
t=text(-0.1,-.1,'$\Gamma^{in}$','Interpreter','latex')
t.FontSize = 14
t=text(-1.9,1.5,'$\Gamma^{out}$','Interpreter','latex')
t.FontSize = 14
print('sqkite_domain.eps','-depsc2')