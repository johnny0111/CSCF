%%
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
%% setup & model
m = 1500;
Cd = 0.3;
area = 1.8;
rho = 1.225;
g = 9.8065;
theta = 0.05;
fmax = 2800;
fmin = -2200;
deltaFmax = 200;%300;
maxTheta = 0.1;
dDist = 10;
distMin = 3;
vRe = 25;
ve = 25;
we = 0;
Ts = 0.1;
d = (rho*area*Cd*ve)/m;

Ac = [0 -1 0 0; 0 -d 0 0; 0 1 0 -1; 0 0 0 -d];
Bc=[0 0; 1/m 0; 0 0; 0 1/m];

A1c = [0 -1;0 -d];
A2c = [-d 0 0;1 0 -1 ; 0 0 -d];
B11c = [0 1/m]';
B12c = [0 0]';
B21c = [1/m 0 0]';
B22c = [0 0 1/m]';
C1 = [1 0];
C2 = [0 1 0];

A = eye(4) + Ac*Ts;
B = Bc*Ts;
C = [1 0 0 0; 0 0 1 0];
A1 = eye(2) + A1c*Ts;
A2 = eye(3) + A2c*Ts;
B11 = B11c*Ts;
B12 = B12c*Ts;
B21 = B21c*Ts;
B22 = B22c*Ts;

x10 = [0;0];
x20 = [0;0;0];
xd0 = [0;0;0;0];

nu = size(B,2);
ny = size(C,1);
nx = size(B,1);
nx1 = size(B11,1);
nu1 = size(B11,2);
nx2 = size(B21,1);
nu2 = size(B21,2);
ny1 = size(C1,1);
ny2 = size(C2,1);


A_ext = [A zeros(nx,ny);C*A eye(ny)];
B_ext = [B;C*B];
C_ext = [zeros(ny,nx),eye(ny)];
A1_ext = [A1 zeros(nx1,ny1);C1*A1 eye(ny1)];
A2_ext = [A2 zeros(nx2,ny2);C2*A2 eye(ny2)];
B11_ext = [B11;C1*B11];
B12_ext = [B12;C1*B12];
B21_ext = [B21;C2*B21];
B22_ext = [B22;C2*B22];
C1_ext = [zeros(ny1,nx1),eye(ny1)];
C2_ext = [zeros(ny2,nx2),eye(ny2)];
%% Cost
% cost parameters
N = 60;
Pi = 15000;
Qi = 0.05;
Ri = 0.001;
alpha1 = 1;
alpha2 = 1;
% distributed steps parameters
w1 = 0.5;
w2 = 1-w1;
np=2;
% Ai={A1_ext,A2_ext};
% Bi={B11_ext,B12_ext;B21_ext,B22_ext};
% Ci={C1_ext,C2_ext};
% [Fb,Gb,Qb,Rb] = GetBatchYNdistMatrices(Ai,Bi,Ci,N,Pi,Qi,Ri,alpha1)
%% Distributed Batch Matrices

[Fb1,Gb11,Gb12,Qb1,Rb1] = GetBatchYdistMatrices(A1_ext,B11_ext,B12_ext,C1_ext,N,Pi,Qi,Ri,alpha1);
[Fb2,Gb22,Gb21,Qb2,Rb2] = GetBatchYdistMatrices(A2_ext,B22_ext,B21_ext,C2_ext,N,Pi,Qi,Ri,alpha2);
[Rt1,S11x,S11y,S12x,S12y,S12u,Rt2,S22x,S22y,S21x,S21y,S21u] = GetDistMPCbatch(Fb1,Gb11,Gb12,Qb1,Rb1,Fb2,Gb21,Gb22,Qb2,Rb2);


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
My1 = [-Gb11; Gb11];

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
My2 = [-Gb22; Gb22];

M2 = [Mdu2; Mu2];






%% Simulation


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
%Y1(:,1) = C1*xd0(1:2);
%Xd1(:,2) = xd0(1:2);
%X1(:,2) = x01;
%Y1(:,2) = C1*xd0(1:2);

%player2
x02 = [xd0(2:4)*0 ; C2*xd0(2:4)];
U2 = zeros(nu2,N,nk);
U2_d = zeros(nu2,N,nk);
dU2 = zeros(nu2,N);
dU2_d = zeros(nu2,N);
%Xd2(:,1) =  xd0(3:5);
X2(:,1) = x02;
Y2(:,1) = C2*xd0(2:4);
%Xd2(:,2) = xd0(3:5);
X2(:,2) = x02;
Y2(:,2) = C2*xd0(2:4);

Xd(:,1) = xd0;
Xd(:,2) = xd0;
u(:,1) = [0;0];
for k = 2:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    
    % distributed MPC
    Dxdk1 = Xd(1:2,k)-Xd(1:2,k-1);
    X1(:,k) = [ Dxdk1; C1*Xd(1:2,k)];
    x1k = X1(:,k);
    
    Dxdk2 = [Xd(2,k);Xd(3,k);Xd(4,k)]-[Xd(2,k-1);Xd(3,k-1);Xd(4,k-1)];
    X2(:,k) = [ Dxdk2; C2*[Xd(2,k);Xd(3,k);Xd(4,k)]];
    x2k = X2(:,k);
    wu1 = [-U_min1 + M41*u(1,k-1);U_max1 - M41*u(1,k-1)];
    wu2 = [-U_min2 + M42*u(2,k-1);U_max2 - M41*u(2,k-1)];
    wy1 = [-Y_min1 + Fb1*x1k; Y_max1 - Fb1*x1k];
    wy2 = [-Y_min2 + Fb2*x2k; Y_max2 - Fb2*x2k];
    
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);

    %wu1 = [-U_min1 + M41*U1p;U_max1 - M41*U1p];
    %wu2 = [-U_min2 + M42*U2p;U_max2 - M42*U2p];
    wr1 = [wdu1;wu1];
    wr2 = [wdu2;wu2];
    % Get optimal sequence for player 1
    St1 = S11x*x1k - S11y*Yb1 + S12x*x2k - S12y*Yb2 + S12u*U2p;
    [U1o,J1o,exitflag,output,lambda] = quadprog(Rt1,St1,M1,wr1);
    if exitflag<0
        error('Problems in the Optimization problem (player 1).');
    end
    
    % Get optimal sequence for player 2
    St2 = S21x*x1k - S21y*Yb1 + S22x*x2k - S22y*Yb2 + S21u*U1p ;
    [U2o,J2o,exitflag,output,lambda] = quadprog(Rt2,St2,M2,wr2);
    if exitflag<0
        error('Problems in the Optimization problem (player 2).');
    end
    
    
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
    for p = 1:np
        U1pp = w1*U1o + w2*U1pp ;
        U2pp = w2*U2o + w1*U2pp ;
    end
    U1(:,:,k) = reshape( U1pp ,nu1,N);
    U2(:,:,k) = reshape( U2pp ,nu2,N);
    % apply first value at each player
    u1k = U1(:,1,k);
    u2k = U2(:,1,k);
    u(:,k) = [ u1k+  u(1,k-1) ; u2k+u(2,k-1) ];
    
    % simulate system for distributed MPC
    Xd(:,k+1) = A*Xd(:,k) + B*u(:,k); %simulate joint system
end


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
legend('car1 cent.','car2 cent');
title('Phase plot');

figure(7102);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
%plot(Tref,ref(2,:),'k+--');
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstdarkblue);
plot(TX,Xd(3,:),'.--','Color',sstlightblue);
plot(TX,Xd(4,:),'-*','Color','magenta');
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','car1 $$x_1$$ cen.','car1 $$x_2$$ cen.','car 2$$x_1$$ cen.','car 2 $$x_2$$ cen.','car 1$$x_1$$ decen.','car 1$$x_2$$  deccen.','car 2$$x_1$$ decen.','car 2$$x_2$$  deccen.','Location','SouthEast');
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
legend('$$u_1$$ cent.','$$u_2$$ cent.', 'du1', 'du2');
title('Input');






























