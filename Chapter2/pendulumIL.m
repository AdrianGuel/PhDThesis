%pendulum stochastic differential equation under the laplace assumption
%Adrian Guel 2022
%Reference of the Laplace assumption "Population dynamics under the laplace
%assumption", Andre C. Marreiros et al. 2009

clearvars;
close all;
clc;

fig=figure('visible','on');
set(fig, 'Position',  [615,328,800,354])
set(gcf,'color','w');
ax1 = subplot(3,2,[2 4 6]);
hold(ax1,'on')
grid(ax1,'on')
xlabel(ax1,'$\theta$','Interpreter','Latex','FontSize', 14)
ylabel(ax1,'$\dot{\theta}$','Interpreter','Latex','FontSize', 14)
axis(ax1,'square')
%view(ax1,3)

ax2 = subplot(3,2,1);
hold(ax2,'on')
grid(ax2,'on')
xlabel(ax2,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax2,'$\Sigma_{11}(t)$','Interpreter','Latex','FontSize', 14)

ax3 = subplot(3,2,3);
hold(ax3,'on')
grid(ax3,'on')
xlabel(ax3,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax3,'$\Sigma_{22}(t)$','Interpreter','Latex','FontSize', 14)

ax4 = subplot(3,2,5);
hold(ax4,'on')
grid(ax4,'on')
xlabel(ax4,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(ax4,'$\Sigma_{12}(t)$','Interpreter','Latex','FontSize', 14)

y0=[5;1;1e-2;0;1e-2];
D=[0,0;0,1e-2];
g=9.81;
L=1;
m=1;
b=1;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
tspan=[0 10];
[t,L] = ode45(@(t,y) LaplacianA(t,y,g,L,b,m,D), tspan, y0, opts);


plot(ax2,t,L(:,3),'k')
plot(ax3,t,L(:,5),'k')
plot(ax4,t,L(:,4),'k')

M=zeros(2,length(t));
G=zeros(1,length(t));
    for k=1:length(t)
        M(:,k)=mvnrnd([L(k,1) L(k,2)],[L(k,3),L(k,4);L(k,4),L(k,5)]);
    end
plot(ax1,M(1,:),M(2,:),'k')


function dydt=LaplacianA(t,y,g,L,b,m,D)
    dydt = zeros(5,1);
    dydt(1)=y(2);
    dydt(2)=-((g*sin(y(1)))/L)+(g*y(3)*sin(y(1)))/(2*L)-b*y(2)/m;
    dydt(3)=2*D(1,1) + 2*y(4);
    dydt(4)=2*D(1,2) + y(5)-(g*y(3)*cos(y(1)))/L-b*y(4)/m;
    dydt(5)=2*D(2,2) - 2*(g*y(4)*cos(y(1)))/L-2*b*y(5)/m;
end
