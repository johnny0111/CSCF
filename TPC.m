A = [1 1;0 1];
B = [0;1];
C=[1 0];
N = 3;
P = eye(2);
Q = eye(2);
R = 4;
x0 = [-4.5;2];
nx = size(B,1);
nu = size(B,2);

% compute the batch matrices
[F,G,Qb,Rb,H] = GetBatchXMatrices(A,B,C,N,P,Q,R);
Qt = F'*Qb*F,
Rt = G'*Qb*G + Rb,
St = G'*Qb*F,