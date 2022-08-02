%Stochastic simulation and histogram computation at a given instant of time
%Guel-Cortez 2022
close all;
clearvars;
clc

%Stochastic simulation
theta=3;
mu=0;
sigma=0.1;
D=(sigma);
x0=1;
N=1e2;
tmax=2;
g = gpuDevice(1);
[t,x]=ornstein_uhlenbeck_euler_maruyama( theta, mu, sigma, x0, tmax, N);

fig = figure;
set(gcf,'color','w');
ax1 = axes();
ax2 = axes('Position',[0.65 0.65 0.28 0.28]);
plot(ax1,t,x)
histogram(ax2,x(:,t==0.4))
grid on
xlabel(ax1,'$t$','Interpreter','Latex','FontSize', 16)
ylabel(ax1,'$x$','Interpreter','Latex','FontSize', 16)
xlabel(ax2,'$x$','Interpreter','Latex','FontSize', 16)
ylabel(ax2,'$f$','Interpreter','Latex','FontSize', 16)
