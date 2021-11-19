Ad = [1 1;0 1];
Bd = [0;1];
Cd=[1 0];
N = 3;
P = 1;
Q = 1;
R = 10;
x0 = [-4.5;2];
nx = size(Bd,1);
nu = size(Bd,2);

A = [Ad zeros(2,1);Cd*Ad 1];
B = [Bd;Cd*Bd];
C = [0 0 1];
% compute the batch matrices
% F=[C;C*A;C*A^2;C*A^3]
% H=diag(C)
[F,G,Qb,Rb,H,Fd,Gd,Hd,Ad,Bd,CD] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
F_ = H*F
G_ = H*G
Q_ = zeros(1,N+1);
R_ = zeros(1,N);

for i = 1:N
    Q_(i) = Q;
end
Q_(N+1) = P;
Q_ = diag(Q_)
for i = 1:N
    R_(i) = R;
end
R_ = diag(R_)
RTil = G_'*Q_*G_ + R_
STil = G_'*Q_

