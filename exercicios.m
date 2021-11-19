clc, clear all;
Ad = [1 1;0 1];
Bd = [0;1];
Cd=[1 0];
N = 3;
P = 1;
Q = 1;
R = 4;
Yb = [0 0 0 0]';
xd0 = [-4.5;2];

[F,G,Qb,Rb,H,Fd,Gd,Hd,A,B,C] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R)

Fb = H*F;
Gb = H*G;
Rt = Gb'*Qb*Gb+Rb;
St = Gb'*Qb;
Ky = Rt^-1 * St;
K = Ky*Fb;

x0 = [xd0-xd0; Cd*xd0]
%X = (F-G*K)*x0 + G*Ky*Yb
dU = -K*x0

xd1 = Ad*xd0 + Bd*dU(1)


x1 = [ xd1-xd0; Cd*xd1]
dU1 = -K*x1 
U = dU(1) + dU1(1)