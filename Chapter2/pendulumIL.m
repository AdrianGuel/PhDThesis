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

fig2=figure('visible','on');
set(fig2, 'Position',  [615,328,800,354])
set(gcf,'color','w');
axa = subplot(1,2,1);
hold(axa,'on')
grid(axa,'on')
xlabel(axa,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(axa,'$\Gamma(t)$','Interpreter','Latex','FontSize', 14)
axis(axa,'square')
%view(ax1,3)

axb = subplot(1,2,2);
hold(axb,'on')
grid(axb,'on')
xlabel(axb,'$t$','Interpreter','Latex','FontSize', 14)
ylabel(axb,'$\mathcal{L}(t)$','Interpreter','Latex','FontSize', 14)

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

dm1=gradient(L(:,1),t);
dm2=gradient(L(:,2),t);
ds11=gradient(L(:,3),t);
ds12=gradient(L(:,4),t);
ds22=gradient(L(:,5),t);
s11=L(:,3);
s12=L(:,4);
s22=L(:,5);
G=Gamma(dm1,dm2,ds11,ds12,ds22,s11,s12,s22);
plot(axa,t,G,'k')
plot(axb,t,cumtrapz(t,G),'k')
function dydt=LaplacianA(t,y,g,L,b,m,D)
    dydt = zeros(5,1);
    dydt(1)=y(2);
    dydt(2)=-((g*sin(y(1)))/L)+(g*y(3)*sin(y(1)))/(2*L)-b*y(2)/m;
    dydt(3)=2*D(1,1) + 2*y(4);
    dydt(4)=2*D(1,2) + y(5)-(g*y(3)*cos(y(1)))/L-b*y(4)/m;
    dydt(5)=2*D(2,2) - 2*(g*y(4)*cos(y(1)))/L-2*b*y(5)/m;
end

function val=Gamma(dm1,dm2,ds11,ds12,ds22,s11,s12,s22)
val=(1./(2*(s12.^2-s11.*s22).^2)).*(ds22.^2.*s11.^2+2.*ds22.*s12.*(-2*ds12.*s11+ds11.*s12) ...
    +2*s12.^2.*(ds12.^2+dm2.*(-dm2.*s11+2.*dm1.*s12))+2.*(s11.*(ds12.^2+dm2.^2.*s11) ...
    -2*(ds11.*ds12+dm1.*dm2.*s11).*s12-dm1.^2.*s12.^2).*s22+(ds11.^2+2*dm1.^2.*s11).*s22.^2);
end
