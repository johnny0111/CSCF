%% Exercise book solutions (CPCS)
%  Module 4: Distributed MPC
%
% Bruno Guerreiro (bj.guerreiro@fct.unl.pt)

clear all;
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


%% Exercise 4.1: basic centralized vs decentralized tracking MPC

% define parameters
A1 = 0.7;
A2 = 0.8;
A = blkdiag(A1,A2);
B11 =  1.1;
B12 =  0.9;
B21 = -0.9;
B22 =  1.1;
B = [B11,B12;B21,B22];
C1 = 1;
C2 = 1;
C = blkdiag(C1,C2);
x10 = 1;
x20 = 1;
x0 = [x10;x20];
nu = size(B,2);
ny = size(C,1);
nx = size(B,1);
nx1 = size(B11,1);
nx2 = size(B22,1);
nu1 = size(B11,2);
nu2 = size(B22,2);

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
[Fb,Gb,Qb,Rb,F,G,H] = GetBatchYMatrices(A,B,C,N,P,Q,R);
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^-1*St;
K = -Ky*Fb;

% compute decentralized tracking controllers
[Fb1,Gb1,Qb1,Rb1,F1,G1,H1] = GetBatchYMatrices(A1,B11,C1,N,Pi,Qi,Ri);
Rt1 = Gb1'*Qb1*Gb1 + Rb1;
St1 = Gb1'*Qb1;
Ky1 = Rt1^(-1)*St1;
K1 = -Ky1*Fb1;
[Fb2,Gb2,Qb2,Rb2,F2,G2,H2] = GetBatchYMatrices(A2,B22,C2,N,Pi,Qi,Ri);
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
U = zeros(nu,nk);
Ud = zeros(nu,nk);
X = zeros(nx,nk+1);
X(:,1) = x0;
Xd = zeros(nx,nk+1);
Xd(:,1) = x0;
for k = 1:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
    % centralized MPC control
    Uopt(:,:,k) = reshape( K*X(:,k)+Ky*Yb ,nu,N);
    Xopt(:,:,k) = reshape( F*X(:,k)+G*(K*X(:,k)+Ky*Yb) ,nx,N+1);
    U(:,k) = Uopt(:,1,k);
    
    % joint system simulation for centralized MPC
    X(:,k+1) = A*X(:,k) + B*U(:,k);
    
    % decentralized MPC control
    U1opt(:,:,k) = reshape( K1*Xd(1,k)+Ky1*Yb1 ,nu1,N);
    X1opt(:,:,k) = reshape( F1*Xd(1,k)+G1*(K1*Xd(1,k)+Ky1*Yb1) ,nx1,N+1);
    U2opt(:,:,k) = reshape( K2*Xd(2,k)+Ky2*Yb2 ,nu2,N);
    X2opt(:,:,k) = reshape( F2*Xd(2,k)+G2*(K2*Xd(2,k)+Ky2*Yb2) ,nx2,N+1);
    Ud(:,k) = [U1opt(:,1,k);U2opt(:,1,k)];   
    
    % joint system simulation for decentralized MPC 
    Xd(:,k+1) = A*Xd(:,k) + B*Ud(:,k);
    
end

figure(7101);
plot(ref(1,:),ref(2,:),'k+-');
grid on;
hold on;
plot(X(1,:),X(2,:),'s-','Color',sstblue);
plot(Xd(1,:),Xd(2,:),'o-','Color',sstgreen);
for k = 1:nk
    plot(Xopt(1,:,k),Xopt(2,:,k),'.--','Color',sstlightblue);
    plot(X1opt(1,:,k),X2opt(1,:,k),'.-.','Color',sstlightgray);
end
plot(X(1,:),X(2,:),'s-','Color',sstblue);
plot(Xd(1,:),Xd(2,:),'o-','Color',sstgreen);
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('ref.','cent.','decent.','pred. cent.','pred. dec.');
title('Phase plot');

figure(7102);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(Tref,ref(2,:),'k+--');
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstdarkblue);
plot(TX,Xd(1,:),'o-','Color',sstgreen);
plot(TX,Xd(2,:),'+-','Color',sstdarkgreen);
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-.','Color',sstlightgray);
end
plot(Tref,ref(1,:),'k+-');
plot(Tref,ref(2,:),'k+--');
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstdarkblue);
plot(TX,Xd(1,:),'o-','Color',sstgreen);
plot(TX,Xd(2,:),'+-','Color',sstdarkgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','$$x_1$$ cen.','$$x_2$$ cen.','$$x_1$$ decen.','$$x_2$$ decen.','$$x_1$$ pred. cen.','$$x_2$$ pred. cen.','Location','SouthEast');
title('State evolution');

figure(7103);
plot(TU,U(1,:),'s-','Color',sstblue);
grid on;
hold on;
plot(TU,U(2,:),'d-','Color',sstdarkblue);
plot(TU,Ud(1,:),'o-','Color',sstgreen);
plot(TU,Ud(2,:),'+-','Color',sstdarkgreen);
for k = 1:nk
    plot(TUopt(k,:),Uopt(1,:,k),'.--','Color',sstlightblue);
    plot(TUopt(k,:),Uopt(2,:,k),'.-.','Color',sstlightgray);
end
plot(TU,U(1,:),'s-','Color',sstblue);
plot(TU,U(2,:),'d-','Color',sstdarkblue);
plot(TU,Ud(1,:),'o-','Color',sstgreen);
plot(TU,Ud(2,:),'+-','Color',sstdarkgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.','$$u_1$$ decent.','$$u_2$$ decent.');
title('Input');

%% Exercise 4.2: no solution provided (centralized vs decentralized tracking incremental MPC)

%% Exercise 4.3: basic distributed tracking MPC

% printFigs = 1;
example_name = 'exercise43_distributed_doubleint_p2';
np = 2;
% example_name = 'exercise43_distributed_doubleint_p10';
% np = 10;

Ai = [  1 , 0.1 , -0.1 
        0 ,   1 ,  0
		0 ,   0 ,  1 ] ;
A1 = Ai;
A2 = Ai;
A = Ai;
B11 = [ 0 ; 1 ; 0];
B22 = [ 0 ; 0 ; 1];
B12 = B22;
B21 = B11;
B = [B11,B22];
C1 = [ 1 0 0];
C2 = [-1 0 0];
C = [C1;C2];
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
N = 3;
Pi = 2;
Qi = 2;
Ri = 1;
alpha1 = 1;
alpha2 = 1;
% distributed steps parameters
w1 = 0.4;
w2 = 1-w1;

% compute centralized tracking controller
P = blkdiag(Pi,Pi);
Q = blkdiag(Qi,Qi);
R = blkdiag(Ri,Ri);
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
% E = blkdiag(E1,E2); % Joint MPC policy selection matrix
Ec = [E1;E2];

% compute distributed tracking controllers
[Fb1,Gb11,Gb12,Qb1,Rb1] = GetBatchYdistMatrices(A1,B11,B12,C1,N,Pi,Qi,Ri,alpha1);
[Fb2,Gb22,Gb21,Qb2,Rb2] = GetBatchYdistMatrices(A2,B22,B21,C2,N,Pi,Qi,Ri,alpha2);
[K11,K12,K21,K22,K11y,K12y,K21y,K22y,L1,L2] = GetDistMPCGains(Fb1,Gb11,Gb12,Qb1,Rb1,Fb2,Gb21,Gb22,Qb2,Rb2);

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
    
    % centralized MPC
    Ucopt(:,:,k) = reshape( Kc*xc(:,k)+Kcy*Yb ,nu,N);
    uc(:,k) = Ucopt(:,1,k);
    
    % simulate system for decentralized MPC
    xc(:,k+1) = A*xc(:,k) + B*uc(:,k); %simulate joint system
    
    % decentralized MPC
    xd1 = xd(:,k);
    xd2 = xd(:,k);
    Ud1opt(:,:,k) = reshape( Kd1*xd1+Kdy1*Yb1 ,nu1,N);
    Ud2opt(:,:,k) = reshape( Kd2*xd2+Kdy2*Yb2 ,nu2,N);
    ud1 = Ud1opt(:,1,k);    
    ud2 = Ud2opt(:,1,k);    
    ud(:,k) = [ ud1 ; ud2 ];
    
    % simulate system for decentralized MPC
    xd(:,k+1) = A*xd(:,k) + B*ud(:,k); %simulate joint system
    
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

figure(7201);
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

figure(7202);
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

figure(7203);
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

%% Exercise 4.4: no solution provided (basic distributed tracking MPC)

%% Exercise 4.5: constrained distributed tracking MPC for 3 agents

clear TU TX Tref ref uc ud u U1 U2 U3 xc xd x Ucopt Ud1opt Ud2opt Ud3opt;
example_name = 'exercise45_constrainedDMPC_3doubleint_p2';

A  = [  1 , 0.1 , -0.1 ,    0 ,   0
        0 ,   1 ,    0 ,    0 ,   0
		0 ,   0 ,    1 ,    0 ,   0 
		0 ,   0 ,    0 ,    1 ,   0 
		0 ,   0 ,  0.1 , -0.1 ,   1 ] ;
Ai = {A(1:3,1:3),A(1:3,1:3),A(3:5,3:5)};
B = [   0 , 0 , 0
        1 , 0 , 0
        0 , 1 , 0 
        0 , 0 , 1 
        0 , 0 , 0 ];
Bij = { B(1:3,1), B(1:3,2), B(1:3,3)
        B(1:3,1), B(1:3,2), B(1:3,3)
        B(3:5,1), B(3:5,2), B(3:5,3) };
C = [   1 , 0 , 0 , 0 , 0
       -1 , 0 , 0 , 0 , 0
        0 , 0 , 0 , 0 , 1 ];
Ci = {C(1,1:3),C(2,1:3),C(3,3:5)};
x0 = [  1 ; 1 ; -1 ; 2 ; 2 ];
x10 = x0(1:3);
x20 = x0(1:3);
x30 = x0(3:5);
nu = size(B,2);
ny = size(C,1);
nx = size(B,1);
nxi = size(Bij{1,1},1);
nui = size(Bij{1,1},2);

% cost parameters
N = 7;
Pi = 1;
Qi = 1;
Ri = 1;
alphai = 1;
% distributed steps parameters
w1 = 0.4;
w2 = 0.3;
w3 = 1-w1-w2;
np = 2;
[Fb,Gb,Qb,Rb] = GetBatchYNdistMatrices(Ai,Bij,Ci,N,Pi,Qi,Ri,alphai);
% compute centralized tracking controller
P = blkdiag(Pi,Pi,Pi);
Q = blkdiag(Qi,Qi,Qi);
R = blkdiag(Ri,Ri,Ri);
[Fbc,Gbc,Qbc,Rbc] = GetYMats(A,B,C,N,P,Q,R);
[Rtc,Sxc,Syc] = GetYMPC(Fbc,Gbc,Qbc,Rbc);

% compute decentralized tracking controllers
[Fbd1,Gbd1,Qbd1,Rbd1] = GetYMats(Ai{1},Bij{1,1},Ci{1},N,Pi,Qi,Ri);
[Rtd1,Sxd1,Syd1] = GetYMPC(Fbd1,Gbd1,Qbd1,Rbd1);
[Fbd2,Gbd2,Qbd2,Rbd2] = GetYMats(Ai{2},Bij{2,2},Ci{2},N,Pi,Qi,Ri);
[Rtd2,Sxd2,Syd2] = GetYMPC(Fbd2,Gbd2,Qbd2,Rbd2);
[Fbd3,Gbd3,Qbd3,Rbd3] = GetYMats(Ai{3},Bij{3,3},Ci{3},N,Pi,Qi,Ri);
[Rtd3,Sxd3,Syd3] = GetYMPC(Fbd3,Gbd3,Qbd3,Rbd3);

% compute distributed tracking controllers
[Fb,Gb,Qb,Rb] = GetYMats(Ai,Bij,Ci,N,Pi,Qi,Ri,alphai);
[Rt,Sx,Sy,Su] = GetYMPC(Fb,Gb,Qb,Rb);

% compute constraints matrices:
Ui_max = 0.5*ones(N*nui,1);
Mui = [-eye(N*nui);eye(N*nui)];
wui = [Ui_max;Ui_max];
U_max = 0.5*ones(N*nu,1);
Mu = [-eye(N*nu);eye(N*nu)];
wu = [U_max;U_max];

%simulate controlled system:
nk = 100;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = [2*(Tref>=50);-2*(Tref>=50);2*(Tref>=50)]; % step reference of amplitude 1, starting at k=30
uc = zeros(nu,nk);
ud = zeros(nu,nk);
u = zeros(nu,nk);
U1 = zeros(nui,N,nk);
U2 = zeros(nui,N,nk);
U3 = zeros(nui,N,nk);
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
    Yb3 = ref(3,k:k+N)';
    
    % centralized
    Stc = Sxc*xc(:,k)-Syc*Yb;
    [Uco,Jco,exitflag] = quadprog(Rtc,Stc,Mu,wu);
    if exitflag~=1, error('Problems in centralized.'); end
    Ucopt(:,:,k) = reshape( Uco ,nu,N);
    uc(:,k) = Ucopt(:,1,k);
    xc(:,k+1) = A*xc(:,k) + B*uc(:,k); %simulate joint system
    
    % decentralized
    Std1 = Sxd1*xd(1:3,k)-Syd1*Yb1;
    Std2 = Sxd2*xd(1:3,k)-Syd2*Yb2;
    Std3 = Sxd3*xd(3:5,k)-Syd3*Yb3;
    [Ud1o,Jd1o,ed1] = quadprog(Rtd1,Std1,Mui,wui);
    [Ud2o,Jd2o,ed2] = quadprog(Rtd2,Std2,Mui,wui);
    [Ud3o,Jd3o,ed3] = quadprog(Rtd3,Std3,Mui,wui);
    if any([ed1,ed2,ed3]~=1), error('Problems in decentralized MPC.'); end
    Ud1opt = reshape( Ud1o ,nui,N);
    Ud2opt = reshape( Ud2o ,nui,N);
    Ud3opt = reshape( Ud3o ,nui,N);
    ud(:,k) = [ Ud1opt(:,1) ; Ud2opt(:,1) ; Ud3opt(:,1) ];
    xd(:,k+1) = A*xd(:,k) + B*ud(:,k); %simulate joint system
    
    % distributed MPC
    x1k = x(1:3,k);
    x2k = x(1:3,k);
    x3k = x(3:5,k);
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    U3p = reshape( U3(:,:,k-1) ,[],1);
    % Get optimal sequence for player 1
    St1 = Sx{1,1}*x1k - Sy{1,1}*Yb1 + Sx{1,2}*x2k - Sy{1,2}*Yb2 + Su{1,2}*U2p;
    [U1o,J1o,exitflag] = quadprog(Rt{1},St1,Mui,wui);
    if exitflag~=1, error('Problems in player 1.'); end
    % Get optimal sequence for player 2
    St2 = Sx{2,1}*x1k - Sy{2,1}*Yb1 + Sx{2,2}*x2k - Sy{2,2}*Yb2 + Su{2,1}*U1p;
    [U2o,J2o,exitflag] = quadprog(Rt{2},St2,Mui,wui);
    if exitflag~=1, error('Problems in player 2.'); end
    % Get optimal sequence for player 3
    St3 = Sx{3,2}*x2k - Sy{3,2}*Yb2 + Sx{3,3}*x3k - Sy{3,3}*Yb3 + Su{3,2}*U2p;
    [U3o,J3o,exitflag] = quadprog(Rt{3},St3,Mui,wui);
    if exitflag~=1, error('Problems in player 3.'); end
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
    U3pp = U3p;
    for p = 1:np
        U1pp = w1*U1o + w2*U1pp;
        U2pp = w2*U2o + w1*U2pp;
        U3pp = w2*U3o + w1*U3pp;
    end
    U1(:,:,k) = reshape( U1pp ,nui,N);
    U2(:,:,k) = reshape( U2pp ,nui,N);
    U3(:,:,k) = reshape( U3pp ,nui,N);
    
    u(:,k) = [ U1(:,1,k) ; U2(:,1,k) ; U3(:,1,k) ]; % apply first input to each player
    x(:,k+1) = A*x(:,k) + B*u(:,k); %simulate joint system
end

figure(7502);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(Tref,ref(3,:),'k.--');
plot(TX,xc(1,:),'s-','Color',sstdarkblue);
plot(TX,xc(5,:),'.-.','Color',sstlightblue);
plot(TX,xd(1,:),'d-','Color',sstdarkgray);
plot(TX,xd(5,:),'.-.','Color',sstgray);
plot(TX,x(1,:),'o-','Color',sstdarkgreen);
plot(TX,x(5,:),'.-.','Color',sstlightgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_3$$','$$y_1$$ cen.','$$y_3$$ cen.','$$y_1$$ decen.','$$y_3$$ decen.','$$y_1$$ dist.','$$y_3$$ dist.','Location','SouthEast');
title('Output evolution');

figure(7503);
plot(TU,uc(1,:),'s-','Color',sstdarkblue);
grid on;
hold on;
plot(TU,uc(2,:),'.--','Color',sstblue);
plot(TU,uc(3,:),'.-.','Color',sstlightblue);
plot(TU,ud(1,:),'d-','Color',sstdarkgray);
plot(TU,ud(2,:),'.--','Color',sstgray);
plot(TU,ud(3,:),'.-.','Color',sstlightgray);
plot(TU,u(1,:),'o-','Color',sstdarkgreen);
plot(TU,u(2,:),'.--','Color',sstgreen);
plot(TU,u(3,:),'.-.','Color',sstlightgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend( '$$u_1$$ cent.','$$u_2$$ cent.','$$u_3$$ cent.',...
        '$$u_1$$ decent.','$$u_2$$ decent.','$$u_3$$ decent.',...
        '$$u_1$$ dist.','$$u_2$$ dist.','$$u_3$$ dist.');
title('Input');

%% Exercise 4.6: no solution provided (constrained distributed tracking MPC)

