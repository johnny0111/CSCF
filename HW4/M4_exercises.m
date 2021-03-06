clear all, clc;
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

%% Exercise 4.2: Centralized vs Decentralized Tracking Incremental MPC

% define parameters
A1 = 0.6;
A2 = 0.9;
A = blkdiag(A1,A2);
B11 =  1.2;
B12 =  0.8;
B21 = -0.8;
B22 =  1.2;
B = [B11,B12;B21,B22];
C1 = 1;
C2 = 1;
C = blkdiag(C1,C2);

x10 = 1;
x20 = 1;
xd0 = [x10;x20];
nu = size(B,2);
ny = size(C,1);
nxd = size(B,1);
nx1 = size(B11,1);
nx2 = size(B22,1);
nu1 = size(B11,2);
nu2 = size(B22,2);
nx = nxd+nu;

N = 3;
Pi = 5;
Qi = 5;
Ri = 1;
x10 = 1;
x20 = 1;

% compute centralized tracking controller
P = blkdiag(Pi,Pi);
Q = blkdiag(Qi,Qi);
R = blkdiag(Ri,Ri);
[Fb,Gb,Qb,Rb,F,G,H,Fd,Gd,~,Ai,Bi,Ci] = GetBatchYiMatrices(A,B,C,N,P,Q,R);
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^-1*St;
K = -Ky*Fb;

% compute decentralized tracking controllers
[Fb1,Gb1,Qb1,Rb1,F1,G1,H1,~,~,~,Ai1,Bi1,~] = GetBatchYiMatrices(A1,B11,C1,N,Pi,Qi,Ri);
Rt1 = Gb1'*Qb1*Gb1 + Rb1;
St1 = Gb1'*Qb1;
Ky1 = Rt1^(-1)*St1;
K1 = -Ky1*Fb1;
[Fb2,Gb2,Qb2,Rb2,F2,G2,H2,~,~,~,Ai2,Bi2,~] = GetBatchYiMatrices(A2,B22,C2,N,Pi,Qi,Ri);
Rt2 = Gb2'*Qb2*Gb2 + Rb2;
St2 = Gb2'*Qb2;
Ky2 = Rt2^(-1)*St2;
K2 = -Ky2*Fb2;

%simulate controlled system:
nk = 60;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = [1*(Tref>=30);1.5*(Tref>=30)]; % step reference of amplitude 1, starting at k=30
x0 = [xd0 ; C*xd0];
U = zeros(nu,nk);
Udec = zeros(nu,nk); % solução temporária
X = zeros(nx,nk+1);
X(:,1) = x0;
Xdec = zeros(nx,nk+1);
Xdec(:,1) = x0;
Xdc = zeros(nxd,nk+1);
Xdc(:,1) = xd0;
Xd = zeros(nxd,nk+1);
Xd(:,1) = xd0;
for k = 2:nk
    % compute initial conditions
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    Dxdck = Xdc(:,k)-Xdc(:,k-1);
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdck; C*Xdc(:,k)]; % new x0 - for centralized
    Xdec(:,k) = [ Dxdk; C*Xd(:,k)]; % new x0 - for decentralized
    
    % centralized MPC control
    dUopt(:,:,k) = reshape( K*X(:,k)+Ky*Yb ,nu,N);
    Xopt(:,:,k) = reshape( F*X(:,k)+G*(K*X(:,k)+Ky*Yb) ,nx,N+1); 
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    U(:,k) = Uopt(:,1,k);
%     Xd_opt(:,:,k) = reshape( Fd*Xdc(:,k)+Gd*Uopt(:,:,k)' ,nxd,N+1);
    
    % joint system simulation for centralized MPC
%     X(:,k+1) = Ai*X(:,k) + Bi*U(:,k);
    Xdc(:,k+1) = A*Xdc(:,k) + B*U(:,k);
    
    % decentralized MPC control
    dU1opt(:,:,k) = reshape( K1*[Xdec(1,k);Xdec(3,k)]+Ky1*Yb1 ,nu1,N); % meter nu1 onde está []
    X1opt(:,:,k) = reshape( F1*[Xdec(1,k);Xdec(3,k)]+G1*(K1*[Xdec(1,k);Xdec(3,k)]+Ky1*Yb1) ,nx1,[]); % meter N+1 onde está []
    dU2opt(:,:,k) = reshape( K2*[Xdec(2,k);Xdec(4,k)]+Ky2*Yb2 ,nu2,N); % meter nu2 onde está []
    X2opt(:,:,k) = reshape( F2*[Xdec(2,k);Xdec(4,k)]+G2*(K2*[Xdec(2,k);Xdec(4,k)]+Ky2*Yb2) ,nx2,[]); % meter N+1 onde está []
    U1opt(:,:,k) = Udec(1,k-1) + dU1opt(:,:,k);
    U2opt(:,:,k) = Udec(2,k-1) + dU2opt(:,:,k);
    Udec(:,k) = [U1opt(:,1,k);U2opt(:,1,k)];   
    
    % joint system simulation for decentralized MPC 
%     Xdec(:,k+1) = Ai*Xdec(:,k) + Bi*Udec(:,k);
    Xd(:,k+1) = A*Xd(:,k) + B*Udec(:,k);
    
end

figure(1);
plot(ref(1,:),ref(2,:),'k+-');
grid on;
hold on;
plot(Xdc(1,:),Xdc(2,:),'s-','Color',sstblue);
plot(Xd(1,:),Xd(2,:),'o-','Color',sstgreen);
for k = 1:nk
%     plot(Xopt(1,:,k),Xopt(2,:,k),'.--','Color',sstlightblue); % original
%     plot(Xd_opt(1,:,k),Xd_opt(2,:,k),'.-','Color',sstlightblue); % from me
%     plot(X1opt(1,1:N+1,k),X2opt(1,1:N+1,k),'.-.','Color',sstlightgray);
end
plot(Xdc(1,:),Xdc(2,:),'s-','Color',sstblue);
plot(Xd(1,:),Xd(2,:),'o-','Color',sstgreen);
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('ref.','cent.','decent.','pred. cent.','pred. dec.');
title('Phase plot');

figure(2);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(Tref,ref(2,:),'k+--');
plot(TX,Xdc(1,:),'s-','Color',sstblue);
plot(TX,Xdc(2,:),'d-','Color',sstdarkblue);
plot(TX,Xd(1,:),'o-','Color',sstgreen);
plot(TX,Xd(2,:),'+-','Color',sstdarkgreen);
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-.','Color',sstlightgray);
end
plot(Tref,ref(1,:),'k+-');
plot(Tref,ref(2,:),'k+--');
plot(TX,Xdc(1,:),'s-','Color',sstblue);
plot(TX,Xdc(2,:),'d-','Color',sstdarkblue);
plot(TX,Xd(1,:),'o-','Color',sstgreen);
plot(TX,Xd(2,:),'+-','Color',sstdarkgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','$$x_1$$ cen.','$$x_2$$ cen.','$$x_1$$ decen.','$$x_2$$ decen.','$$x_1$$ pred. cen.','$$x_2$$ pred. cen.','Location','SouthEast');
title('State evolution');

figure(3);
plot(TU,U(1,:),'s-','Color',sstblue);
grid on;
hold on;
plot(TU,U(2,:),'d-','Color',sstdarkblue);
plot(TU,Udec(1,:),'o-','Color',sstgreen);
plot(TU,Udec(2,:),'+-','Color',sstdarkgreen);
for k = 1:nk
    plot(TUopt(k,:),Uopt(1,:,k),'.--','Color',sstlightblue);
    plot(TUopt(k,:),Uopt(2,:,k),'.-.','Color',sstlightgray);
end
plot(TU,U(1,:),'s-','Color',sstblue);
plot(TU,U(2,:),'d-','Color',sstdarkblue);
plot(TU,Udec(1,:),'o-','Color',sstgreen);
plot(TU,Udec(2,:),'+-','Color',sstdarkgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.','$$u_1$$ decent.','$$u_2$$ decent.');
title('Input');

%% Exercise 4.4: Basic Distributed Tracking MPC

np = 2;

% system
A = [0.5   0       0.1    -0.2    0 
      0     0.3      0     -0.1    0.9          
      0     0       0.8    0       0 
      0     0       0       1       0
      0     0       0       0       0.7] ;
A1 = [0.5 0.1 -0.2
      0   0.8  0
      0   0    1];
  
A2 = [0.3   -0.1    0.9
      0     1       0
      0     0       0.7];
A3 = A2;
Ai={A1,A2,A3};

B=[0 0 0;0 0 0;1 0 0; 0 1 0;0 0 1];
B11 = [ 0 ; 1 ; 0];
B12 = [ 0 ; 0 ; 1];
B13 = [0;0;0];
B21 = [0;0;0];
B22=[0;1;0];
B23=[0;0;1];
B31=B21;
B32=B22;
B33=B23;
Bij={B11, B12, B13
    B21,  B22,  B23
    B31,  B32,  B33};
C=[1 0 0 0 0;0 1 0 0 0];
C1 = [ 1 0 0];
C2 = [0 -1 0];
C3 = [0 1 0];
Ci={C1,C2,C3};
x10 = [1;1;2];
x20 = x10;
x0 = x10;
nu = size(B,2);
ny = size(C,1);
nx = size(B,1);
nx1 = size(B11,1);
nx2 = size(B22,1);
nu1 = size(B11,2);
nu2 = size(B22,2);

% cost parameters
N = 2; %3
Pi = 1; %1
Qi = 1; %1
Ri = 1;
alpha1 = 1;
alpha2 = 1;
% distributed steps parameters
w1 = 0.4;
w2 = 1-w1;



% compute distributed tracking controllers
[Fb,Gb,Qb,Rb] = GetBatchYNdistMatrices(Ai,Bij,Ci,N,Pi,Qi,Ri,1);
[Rt,Sx,Sy,Su,K,Ky,Li] = GetYMPC(Fb,Gb,Qb,Rb)
E1 = [eye(nu1) zeros(nu1,nu1*(N-1))]
E2 = [eye(nu2) zeros(nu2,nu2*(N-1))]
% analyze stability of distributed law:
L = [ (1-w1)*eye(size(L1))  , w1*L1
      w2*L2                 , (1-w2)*eye(size(L2)) ];
maxL = max(abs(eig(L))),
% Kb = [K11,K12;K21,K22];
% Kby = [K11y,K12y;K21y,K22y];
% K = E*(1-L)^(-1)*Kb;
% lbd_cl = eig(A+B*K),

%simulate controlled system:
nk = 100;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = [2*(Tref>=50);-2*(Tref>=50)]; % step reference of amplitude 1, starting at k=30
uc = zeros(nu,nk);
ud = zeros(nu,nk);
u = zeros(nu,nk);
U1 = zeros(nu1,N,nk);
U2 = zeros(nu2,N,nk);
xc = zeros(nx,nk+1);
xd = zeros(nx,nk+1);
x = zeros(nx,nk+1);
xc(:,1) = x0;
xd(:,1) = x0;
x(:,1) = x0;
xc(:,2) = x0;
xd(:,2) = x0;
x(:,2) = x0;
for k = 2:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    
    
    % distributed MPC
    x1k = x(:,k);
    x2k = x(:,k);
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    % Get optimal sequence for player 1
    U1o = K11*x1k + K11y*Yb1 + K12*x2k + K12y*Yb2 + L1*U2p;
    % Get optimal sequence for player 2
    U2o = K21*x1k+K21y*Yb1+K22*x2k+K22y*Yb2+L2*U1p;
    % convex step iterates with comms between players
    U1pp = U1p;
    U2pp = U2p;
    for p = 1:np
        U1pp = w1*U1o + (1-w1)*U1pp;
        U2pp = w2*U2o + (1-w2)*U2pp;
    end
    U1(:,:,k) = reshape( U1pp ,nu1,N);
    U2(:,:,k) = reshape( U2pp ,nu2,N);
    % apply first value at each player
    u1k = U1(:,1,k);
    u2k = U2(:,1,k);
    
    % simulate system for distributed MPC
    u(:,k) = [ u1k ; u2k ];
    x(:,k+1) = A*x(:,k) + B*u(:,k); %simulate joint system
end

figure(4);
plot(ref(1,:),ref(2,:),'k+-');
grid on;
hold on;
plot(xc(1,:),xc(2,:),'s-','Color',sstblue);
plot(xc(1,:),xc(3,:),'.-','Color',sstdarkblue);
plot(xd(1,:),xd(2,:),'o-','Color',sstgray);
plot(xd(1,:),xd(3,:),'.-','Color',sstlightgray);
plot(x(1,:),x(2,:),'x-','Color',sstgreen);
plot(x(1,:),x(3,:),'.--','Color',sstdarkgreen);
hold off;
xlabel('$$\theta_{12}$$');
ylabel('$$\omega_i$$');
legend('ref.','cent. 1','cent. 2','decent. 1','decent. 2','dist. 1','dist. 2');
title('Phase plot');

figure(5);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(Tref,-ref(2,:),'k.--');
plot(TX,xc(1,:),'s-','Color',sstdarkblue);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'.-.','Color',sstlightblue);
plot(TX,xd(1,:),'d-','Color',sstdarkgray);
plot(TX,xd(2,:),'.--','Color',sstgray);
plot(TX,xd(3,:),'.-.','Color',sstgray);
plot(TX,x(1,:),'o-','Color',sstdarkgreen);
plot(TX,x(2,:),'.--','Color',sstgreen);
plot(TX,x(3,:),'.-.','Color',sstlightgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','$$x_1$$ cen.','$$x_2$$ cen.','$$x_3$$ cen.','$$x_1$$ decen.','$$x_2$$ decen.','$$x_3$$ decen.','$$x_1$$ dist.','$$x_2$$ dist.','$$x_3$$ dist.');
title('State evolution');

figure(6);
plot(TU,uc(1,:),'s-','Color',sstdarkblue);
grid on;
hold on;
plot(TU,uc(2,:),'.--','Color',sstblue);
plot(TU,ud(1,:),'d-','Color',sstdarkgray);
plot(TU,ud(2,:),'.-.','Color',sstgray);
plot(TU,u(1,:),'o-','Color',sstdarkgreen);
plot(TU,u(2,:),'.--','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.','$$u_1$$ decent.','$$u_2$$ decent.','$$u_1$$ dist.','$$u_2$$ dist.');
title('Input');

%% Exercise 4.6: Constrained Distributed Tracking MPC

% system
Ai = [  0.7 , 0.2 , -0.3 ,   0 ,   0 
          0 ,   1 ,    0 ,   0 ,   0
		  0 ,   0 ,  0.9 ,   0 ,   0
          0 ,   0 ,    0 ,   1 ,   0
          0 ,   0 ,  0.2 , -0.3, 0.6] ;
A1 = Ai;
A2 = Ai;
A3 = Ai;
A = Ai;
B11 = [ 0 ; 1 ; 0; 0; 0];
B22 = [ 0 ; 0 ; 1; 0; 0];
B33 = [ 0 ; 0 ; 0; 1; 0];
B12 = B22; % que fazer aqui????
B13 = B33;
B21 = B11;
B23 = B33;
B31 = B11;
B32 = B22;
B = [B11,B22,B33];
C1 = [ 1 0 0 0 0];
C2 = [0 -1 0 0 0];
C3 = [0  0 0 0 1];
C = [C1;C2;C3];
x10 = [1;1;-1;1;2];
x20 = x10;
x30 = x10;
x0 = x10;
nu = size(B,2);
ny = size(C,1);
nx = size(B,1);
nx1 = size(B11,1);
nx2 = size(B22,1);
nx3 = size(B33,1);
nu1 = size(B11,2);
nu2 = size(B22,2);
nu3 = size(B33,2);

% cost parameters
N = 10; % mudar para 10------------------------------------------------------------------------
Pi = 2;
Qi = 2;
Ri = 0.1;
alpha1 = 1;
alpha2 = 1;
alpha3 = 1;
% distributed steps parameters
w1 = 0.4; % donde vem estes valores?????
w2 = 0.4;
w3 = 1-w2-w1; 
np = 2;

% compute centralized tracking controller
P = blkdiag(Pi,Pi,Pi);
Q = blkdiag(Qi,Qi,Qi);
R = blkdiag(Ri,Ri,Ri);
[Fb,Gb,Qb,Rb,F,G,H] = GetBatchYMatrices(A,B,C,N,P,Q,R);
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Kcy = Rt^-1*St;
Kc = -Kcy*Fb;

% compute decentralized tracking controllers
[Fb1,Gb1,Qb1,Rb1,F1,G1,H1] = GetBatchYMatrices(A1,B11,C1,N,Pi,Qi,Ri);
Rt1 = Gb1'*Qb1*Gb1 + Rb1;
St1 = Gb1'*Qb1;
Kdy1 = Rt1^(-1)*St1;
Kd1 = -Kdy1*Fb1;
E1 = [eye(nu1) zeros(nu1,nu1*(N-1))]; %MPC policy selection matrix
[Fb2,Gb2,Qb2,Rb2,F2,G2,H2] = GetBatchYMatrices(A2,B22,C2,N,Pi,Qi,Ri);
Rt2 = Gb2'*Qb2*Gb2 + Rb2;
St2 = Gb2'*Qb2;
Kdy2 = Rt2^(-1)*St2;
Kd2 = -Kdy2*Fb2;
E2 = [eye(nu2) zeros(nu2,nu2*(N-1))]; %MPC policy selection matrix
[Fb3,Gb3,Qb3,Rb3,F3,G3,H3] = GetBatchYMatrices(A3,B33,C3,N,Pi,Qi,Ri);
Rt3 = Gb3'*Qb3*Gb3 + Rb3;
St3 = Gb3'*Qb3;
Kdy3 = Rt3^(-1)*St3;
Kd3 = -Kdy3*Fb3;
E3 = [eye(nu3) zeros(nu3,nu3*(N-1))]; %MPC policy selection matrix
% E = blkdiag(E1,E2); % Joint MPC policy selection matrix
Ec = [E1;E2;E3];

% compute distributed tracking controllers
[Fb1,Gb11,Gb12,Gb13,Qb1,Rb1] = GetBatchYmdistMatrices(A1,B11,B12,B13,C1,N,Pi,Qi,Ri,alpha1);
[Fb2,Gb22,Gb21,Gb23,Qb2,Rb2] = GetBatchYmdistMatrices(A2,B22,B21,B23,C2,N,Pi,Qi,Ri,alpha2);
[Fb3,Gb33,Gb31,Gb32,Qb3,Rb3] = GetBatchYmdistMatrices(A3,B33,B31,B32,C3,N,Pi,Qi,Ri,alpha3);
[Rt1,S11x,S11y,S12x,S12y,S13x,S13y,S12u,S13u,Rt2,S22x,S22y,S21x,S21y,S23x,S23y,S21u,S23u,Rt3,S33x,S33y,S31x,S31y,S32x,S32y,S31u,S32u] = ...
    GetmDMPCbatch(Fb1,Gb11,Gb12,Gb13,Qb1,Rb1,Fb2,Gb22,Gb21,Gb23,Qb2,Rb2,Fb3,Gb33,Gb31,Gb32,Qb3,Rb3);

% compute constraints matrices:
u_max = 0.5;
u_min = -0.4;
U_max = kron(u_max,ones(N,1));
U_min = kron(u_min,ones(N,1));
Mui = [-eye(N*nu1);eye(N*nu1)];
wui = [-U_min;U_max];
Mi = [Mui];
wi = [wui];

%simulate controlled system:
nk = 150;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = [2*(Tref>=50);-2*(Tref>=50);2*(Tref>=50)]; % step reference of amplitude 1, starting at k=30
uc = zeros(nu,nk);
ud = zeros(nu,nk);
u = zeros(nu,nk);
U1 = zeros(nu1,N,nk);
U2 = zeros(nu2,N,nk);
U3 = zeros(nu3,N,nk);
% tive de acrescentar desde aqui --
Ucopt = zeros(nu,N,nk);
Ud1opt = zeros(nu1,N,nk);
Ud2opt = zeros(nu2,N,nk);
Ud3opt = zeros(nu3,N,nk);
% -- até aqui
xc = zeros(nx,nk+1);
xd = zeros(nx,nk+1);
x = zeros(nx,nk+1);
xc(:,1) = x0;
xd(:,1) = x0;
x(:,1) = x0;
xc(:,2) = x0;
xd(:,2) = x0;
x(:,2) = x0;
xc(:,3) = x0;
xd(:,3) = x0;
x(:,3) = x0;
for k = 2:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    Yb3 = ref(3,k:k+N)';
    
    % centralized
    Ucopt(:,:,k) = reshape( Kc*xc(:,k)+Kcy*Yb ,nu,N);
    uc(:,k) = Ucopt(:,1,k);
    
    % simulate system for centralized MPC
    xc(:,k+1) = A*xc(:,k) + B*uc(:,k); %simulate joint system
    
    % decentralized
    xd1 = xd(:,k);
    xd2 = xd(:,k);
    xd3 = xd(:,k);
    Ud1opt(:,:,k) = reshape( Kd1*xd1+Kdy1*Yb1 ,nu1,N);
    Ud2opt(:,:,k) = reshape( Kd2*xd2+Kdy2*Yb2 ,nu2,N);
    Ud3opt(:,:,k) = reshape( Kd3*xd3+Kdy3*Yb3 ,nu3,N);
    ud1 = Ud1opt(:,1,k);    
    ud2 = Ud2opt(:,1,k);    
    ud3 = Ud3opt(:,1,k);
    ud(:,k) = [ ud1 ; ud2 ; ud3 ];
    
    % simulate system for decentralized MPC
    xd(:,k+1) = A*xd(:,k) + B*ud(:,k); %simulate joint system
    
    % distributed MPC
    x1k = x(:,k);
    x2k = x(:,k);
    x3k = x(:,k);
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    U3p = reshape( U3(:,:,k-1) ,[],1);
    % Get optimal sequence for player 1
    St1 = S11x*x1k - S11y*Yb1 + S12x*x2k - S12y*Yb2 + S13x*x3k - S13y*Yb3 + S12u*U2p + S13u*U3p;
    [U1o,J1o,exitflag,output,lambda] = quadprog(Rt1,St1,Mi,wi);
    if exitflag~=1
        error('Problems in the Optimization problem (player 1).');
    end
    % Get optimal sequence for player 2
    St2 = S21x*x1k - S21y*Yb1 + S22x*x2k - S22y*Yb2 + S23x*x3k - S23y*Yb3 + S21u*U1p + S23u*U3p;
    [U2o,J2o,exitflag,output,lambda] = quadprog(Rt2,St2,Mi,wi);
    if exitflag~=1
        error('Problems in the Optimization problem (player 2).');
    end
    % Get optimal sequence for player 3
    St3 = S31x*x1k - S31y*Yb1 + S32x*x2k - S32y*Yb2 + S33x*x3k - S33y*Yb3 + S31u*U1p + S32u*U2p;
    [U3o,J3o,exitflag,output,lambda] = quadprog(Rt3,St3,Mi,wi);
    if exitflag~=1
        error('Problems in the Optimization problem (player 3).');
    end
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
    U3pp = U3p;
    for p = 1:np
        U1pp = w1*U1o + w2*U1pp + w3*U1pp;
        U2pp = w2*U2o + w1*U2pp + w3*U2pp;
        U3pp = w3*U3o + w1*U3pp + w2*U3pp;
    end
    U1(:,:,k) = reshape( U1pp ,nu1,N);
    U2(:,:,k) = reshape( U2pp ,nu2,N);
    U3(:,:,k) = reshape( U3pp ,nu3,N);
    % apply first value at each player
    u1k = U1(:,1,k);
    u2k = U2(:,1,k);
    u3k = U3(:,1,k);
    u(:,k) = [ u1k ; u2k ; u3k ];
    
    % simulate system for distributed MPC
    x(:,k+1) = A*x(:,k) + B*u(:,k); %simulate joint system
end

figure(7);
plot(ref(1,:),ref(2,:),'k+-');
grid on;
hold on;
plot(xc(1,:),xc(2,:),'s-','Color',sstdarkblue);
plot(xc(1,:),xc(3,:),'.--','Color',sstblue);
plot(xd(1,:),xd(2,:),'d-','Color',sstdarkgray);
plot(xd(1,:),xd(3,:),'.-.','Color',sstgray);
plot(x(1,:),x(2,:),'o-','Color',sstdarkgreen);
plot(x(1,:),x(3,:),'.--','Color',sstgreen);
hold off;
xlabel('$$\theta_{12}$$');
ylabel('$$\omega_i$$');
legend('ref.','cent. 1','cent. 2','decent. 1','decent. 2','dist. 1','dist. 2');
title('Phase plot');

figure(8);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(Tref,-ref(2,:),'k.--');
plot(TX,xc(1,:),'s-','Color',sstdarkblue);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'.-.','Color',sstlightblue);
plot(TX,xd(1,:),'d-','Color',sstdarkgray);
plot(TX,xd(2,:),'.--','Color',sstgray);
plot(TX,xd(3,:),'.-.','Color',sstgray);
plot(TX,x(1,:),'o-','Color',sstdarkgreen);
plot(TX,x(2,:),'.--','Color',sstgreen);
plot(TX,x(3,:),'.-.','Color',sstlightgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','$$x_1$$ cen.','$$x_2$$ cen.','$$x_3$$ cen.','$$x_1$$ decen.','$$x_2$$ decen.','$$x_3$$ decen.','$$x_1$$ dist.','$$x_2$$ dist.','$$x_3$$ dist.');
title('State evolution');

figure(9);
plot(TU,uc(1,:),'s-','Color',sstdarkblue);
grid on;
hold on;
plot(TU,uc(2,:),'.--','Color',sstblue);
plot(TU,ud(1,:),'d-','Color',sstdarkgray);
plot(TU,ud(2,:),'.-.','Color',sstgray);
plot(TU,u(1,:),'o-','Color',sstdarkgreen);
plot(TU,u(2,:),'.--','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.','$$u_1$$ decent.','$$u_2$$ decent.','$$u_1$$ dist.','$$u_2$$ dist.');
title('Input');