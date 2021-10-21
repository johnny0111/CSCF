%%
clear all;
clc;
m = 1500;
Cd = 0.3;
area = 1.8;
rho = 1.225;
g = 9.8065;
theta = 0.05;
fmax = 2600;
fmin = -2000;
deltaFmax = 50;
maxTheta = 0.1;
dist = 10;
distMin = 5;
vRe = 20;
ve = 20;
we = 0;
T = 0.1;

%%
%Model
A = [0 -1; 0 (rho*area*Cd*ve) / m];
B=[0; 1/m];
Bdist = [1 0 0;0 -g -(rho*area*Cd*2*ve)/(2*m)];
C=[1 0];
D = 0;
Ad = eye(2) + A*T;
Bd = B*T;
Bdistd = Bdist * T;
lambda = eig(Ad);
%%
%simulation
sys = ss(Ad,Bd,eye(2),D,T);
x0 = [-0.2 0.3];
t = 0:T:500;
u = 30*ones(length(t),1);
lsim(sys,u,t,x0)
grid on

















