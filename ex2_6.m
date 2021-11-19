clear all, clc;
A = [1 1; 0 1];
B=[0;1];
N = 3;
P = eye(2);
Q = P;
R=10;
x0 = [-4.5;2];

[F,G,Qb,Rb,H] = GetBatchXMatrices(A,B,[],N,P,Q,R)

Rt = G'*Qb*G + Rb;
St = G'*Qb*F;
Qt = F'*Qb*F;

