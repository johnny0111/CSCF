%%
clear all;
close all;
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
dDist = 10;
distMin = 5;
vRe = 20;
ve = 20;
we = 0;
T = 0.1;

%%
A1c = [0 -1; 0 (rho*area*Cd*ve) / m];
B1c=[0; 1/m];
C1=[1 0];
Ad = eye(2) + A1c*T;
Bd = B1c*T;

