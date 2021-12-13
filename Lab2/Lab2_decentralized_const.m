%%
% function flag = Lab2_CD_const(Ri,Qi)
clear all;
close all;
clc;

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstdarkgreen    = [20,120,60]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstdarkgray     = [150,150,150]/255;
sstgray         = [190,190,190]/255;
sstlightgray    = [230,230,230]/255;

m = 1500;
Cd = 0.3;
area = 1.8;
rho = 1.225;
g = 9.8065;
theta = 0.05;
fmax = 2800;
fmin = -2200;
deltaFmax = 200%300;
maxTheta = 0.1;
dDist = 10;
distMin = 3;
vRe = 25;
ve = 25;
we = 0;
Ts = 0.1;

N=75;

% Qi=10;%2600
% Ri=0.0002*Qi;%0.00017;
Qi = 153;%100
Ri = 0.000016;
Pi=100;
% Pi = 10;
Qi =38;%100
Ri = 0.0001%0.00019;
Pi=50000;
N=50;

P = blkdiag(Pi,Pi);
Q = blkdiag(Qi,Qi);
R = blkdiag(Ri,Ri);

nk = 600;%250;
TU = 1:nk;
TX = 1:nk+1;
%TX = 0:T:(nk+1)*T - T;
Tref = 1:nk+N;
ref = 10 * square(0.0002*Tref, 0.79);
ref = [ref; ref];
%%
%decentralized
Ac = [0 -1 0 0; 0 -(rho*area*Cd*ve)/m 0 0; 0 1 0 -1; 0 0 0 -(rho*area*Cd*ve)/m];
Bc=[0 0; 1/m 0; 0 0; 0 1/m];
A1c = [0 -1; 0 -(rho*area*Cd*ve)/m];
A2c = [0 -1; 0 -(rho*area*Cd*ve)/m];
B1c = [0; 1/m];
B2c = [0; 1/m];
C1 = [1 0];
C2 = [1 0];
A1 = eye(2) + A1c*Ts;
A2 = eye(2) + A2c*Ts;
B1 = B1c*Ts;
B2 = B2c*Ts;

A = eye(4) + Ac*Ts;
B = Bc*Ts;
C = [1 0 0 0; 0 0 1 0];

x10 = [1 0]';
x20 = [1 0]';
xd0 = [x10;x20];

nx1 = size(B1,1);
nx2 = size(B2,1);
nu1 = size(B1,2);
nu2 = size(B2,2);
ny1 = size(C1,1);
ny2 = size(C2,1);

A1_ext = [A1 zeros(nx1,ny1);C1*A1 eye(ny1)];
A2_ext = [A2 zeros(nx2,ny2);C2*A2 eye(ny2)];
B11_ext = [B1;C1*B1];
B22_ext = [B2;C2*B2];
C1_ext = [zeros(ny1,nx1),eye(ny1)];
C2_ext = [zeros(ny2,nx2),eye(ny2)];

[F1,G1,Qb1,Rb1,H1,Fd1,Gd1,Hd1] = GetBatchXiMatrices(A1,B1,C1,N,Pi,Qi,Ri);
Gb1 = H1*G1;
Fb1 = H1*F1;
Rt1 = Gb1'*Qb1*Gb1 + Rb1;
St1 = Gb1'*Qb1;
Ky1 = Rt1^(-1)*St1;
K1 = -Ky1*Fb1;

[F2,G2,Qb2,Rb2,H2,Fd2,Gd2,Hd2] = GetBatchXiMatrices(A2,B2,C2,N,Pi,Qi,Ri);
Gb2 = H2*G2;
Fb2 = H2*F2;
Rt2 = Gb2'*Qb2*Gb2 + Rb2;
St2 = Gb2'*Qb2;
Ky2 = Rt2^(-1)*St2;
K2 = -Ky2*Fb2;


%% Constraints 

%player 1

u_max1 = fmax -(0.5*rho*area*Cd*ve^2);
u_min1 = fmin - (0.5*rho*area*Cd*ve^2);
U_max1 = kron(u_max1,ones(N*nu1,1));
U_min1 = kron(u_min1,ones(N*nu1,1));
M31 = kron(tril(ones(N)), eye(nu1));
M41 = kron(ones(N,1), eye(nu1));
Mu1 = [-M31;M31];

du_max1 = deltaFmax;
du_min1 = -deltaFmax;
DU_max1 = kron(du_max1,ones(N*nu1,1));
DU_min1 = kron(du_min1,ones(N*nu1,1));
Mdu1 = [-eye(N*nu1);eye(N*nu1)];
wdu1 = [-DU_min1;DU_max1];

pr_min = -7 ;
pr_max = 93;
Y_min1 = kron(pr_min,ones(nu1*(N+1),1));
Y_max1 = kron(pr_max, ones(nu1*(N+1),1));
My = [-Gb1; Gb1];

M1 = [Mdu1; Mu1];





%player 2

u_max2 = fmax -(0.5*rho*area*Cd*ve^2);
u_min2 = fmin - (0.5*rho*area*Cd*ve^2);
U_max2 = kron(u_max2,ones(N*nu2,1));
U_min2 = kron(u_min2,ones(N*nu2,1));
M32 = kron(tril(ones(N)), eye(nu2));
M42 = kron(ones(N,1), eye(nu2));
Mu2 = [-M32;M32];

du_max2 = deltaFmax;
du_min2 = -deltaFmax;
DU_max2 = kron(du_max2,ones(N*nu2,1));
DU_min2 = kron(du_min2,ones(N*nu2,1));
Mdu2 = [-eye(N*nu2);eye(N*nu2)];
wdu2 = [-DU_min2;DU_max2];

pr_min = -7 ;
pr_max = 93;
Y_min2 = kron(pr_min,ones(nu2*(N+1),1));
Y_max2 = kron(pr_max, ones(nu2*(N+1),1));
My2 = [-Gb2; Gb2];

M2 = [Mdu2; Mu1];

%% simulation

nk = 600;%250;
TU = 1:nk;
TX = 1:nk+1;
%TX = 0:T:(nk+1)*T - T;
Tref = 1:nk+N;
ref = 10 * square(0.0002*Tref, 0.79);
ref = [ref; ref];

%player1
x01 = [xd0(1:2)*0 ; C1*xd0(1:2)];
U1 = zeros(nu1,N,nk);
U1_d = zeros(nu1,N,nk);
dU1 = zeros(nu1,N);
dU1_d = zeros(nu1,N);
%Xd1(:,1) = xd0(1:2);
X1(:,1) = x01;
Y1(:,1) = C1*xd0(1:2);
Xd1(:,2) = xd0(1:2);
X1(:,2) = x01;
Y1(:,2) = C1*xd0(1:2);

%player2
x02 = [xd0(3:4)*0 ; C2*xd0(3:4)];
U2 = zeros(nu2,N,nk);
U2_d = zeros(nu2,N,nk);
dU2 = zeros(nu2,N);
dU2_d = zeros(nu2,N);
%Xd2(:,1) =  xd0(3:5);
X2(:,1) = x02;
Y2(:,1) = C2*xd0(3:4);
%Xd2(:,2) = xd0(3:5);
X2(:,2) = x02;
Y2(:,2) = C2*xd0(3:4);

Xd(:,1) = xd0;
Xd(:,2) = xd0;
u(:,1) = [0;0];

for k = 2:nk
    
    %?
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';

    
    
    Dxdk1 = Xd(1:2,k)-Xd(1:2,k-1);
    X1(:,k) = [ Dxdk1; C1*Xd(1:2,k)];
    x1k = X1(:,k); 
    
    
    Dxdk2 = Xd(3:4,k)-Xd(3:4,k-1);
    X2(:,k) = [ Dxdk2; C2*Xd(3:4,k)];
    x2k = X2(:,k);
    
    wu1 = [-U_min1 + M41*u(1,k-1);U_max1 - M41*u(1,k-1)];
    wu2 = [-U_min2 + M42*u(2,k-1);U_max2 - M41*u(2,k-1)];
    
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    
    wr1 = [wdu1;wu1];
    wr2 = [wdu2;wu2];
    
    % Decentralized MPC
      [dUo1,Jo1,exitflag,output1,lambda1] = quadprog_v2(2*Rt1,2*St1*(Fb1*x1k-Yb1),M1,wr1);
    if exitflag<0
            k
            error('Problems in the Optimization problem (1)');          

    end    
    
          [dUo2,Jo2,exitflag,output2,lambda2] = quadprog_v2(2*Rt2,2*St2*(Fb2*x2k-Yb2),M2,wr2);
    if exitflag<0
            k
            error('Problems in the Optimization problem (2)');         

    end 
      

    % apply first value at each player
    u1k = dUo1(1);
    u2k = dUo2(1);
    u(:,k) = [ u1k+  u(1,k-1) ; u2k+u(2,k-1) ];
    
    % simulate system for distributed MPC
    Xd(:,k+1) = A*Xd(:,k) + B*u(:,k); %simulate joint system
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];
MAX = u_max1*ones(1,nk);

%%

figure(7101);

grid on;
hold on;
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
plot(Xd(3,:), Xd(4,:), 's-','Color','magenta');

hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('car1 decent.','car2 decent');
title('Phase plot');

figure(7102);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(Tref,ref(2,:),'k+--');
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstdarkblue);
plot(TX,Xd(3,:),'.--','Color',sstlightblue);
plot(TX,Xd(4,:),'-*','Color','magenta');
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','car1 $$x_1$$ decent.','car1 $$x_2$$ decent.','car2$$x_1$$ decent.','car2 $$x_2$$ decent.','car 1$$x_1$$ decen.','car 1$$x_2$$  deccen.','car 2$$x_1$$ decen.','car 2$$x_2$$  deccen.','Location','SouthEast');
title('State evolution');

figure(7103);
plot(TU,u(1,:),'s-','Color','magenta');
grid on;
hold on;
plot(TU,u(2,:),'d-','Color',sstdarkblue);
plot(TU, MAX);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ decent.','$$u_2$$ decent.', 'du1', 'du2');
title('Input');





%end