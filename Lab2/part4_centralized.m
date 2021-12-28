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

Ac = [0  -1  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
     0  -d  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
     0  1   0   -1  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
     0  0   0   -d  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
     0  0   0   1   0   -1  0   0   0   0   0   0   0   0   0   0   0   0   0   0
     0  0   0   0   0   -d  0   0   0   0   0   0   0   0   0   0   0   0   0   0
     0  0   0   0   0   1   0   -1  0   0   0   0   0   0   0   0   0   0   0   0
     0  0   0   0   0   0   0   -d  0   0   0   0   0   0   0   0   0   0   0   0
     0  0   0   0   0   0   0   1   0   -1  0   0   0   0   0   0   0   0   0   0
     0  0   0   0   0   0   0   0   0   -d  0   0   0   0   0   0   0   0   0   0
     0  0   0   0   0   0   0   0   0   1   0   -1  0   0   0   0   0   0   0   0
     0  0   0   0   0   0   0   0   0   0   0   -d  0   0   0   0   0   0   0   0
     0  0   0   0   0   0   0   0   0   0   0   1   0   -1  0   0   0   0   0   0
     0  0   0   0   0   0   0   0   0   0   0   0   0   -d  0   0   0   0   0   0
     0  0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1  0   0   0   0
     0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   -d  0   0   0   0
     0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1  0   0
     0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -d  0   0
     0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1
     0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -d];
%  Ac = [0  -1  0   0   0   0   0    0  0    0    0   0
%        0  -d  0   0   0   0   0    0  0    0    0   0
%        0  1   0   -1  0   0   0    0  0    0    0   0
%        0  0   0   -d  0   0   0    0  0    0    0   0 
%        0  0   0    1  0   -1  0    0  0    0    0   0
%        0  0   0    0  0   -d  0    0  0    0    0   0  
%        0  0   0   0   0   1   0   -1  0    0    0   0
%        0  0   0   0   0   0   0   -d  0    0    0   0
%        0  0   0   0   0   0   0   1   0   -1    0   0   
%        0  0   0   0   0   0   0   0   0   -d    0   0   
%        0  0   0   0   0   0   0   0   0    1     0   -1  
%        0  0   0   0   0   0   0   0   0    0     0   -d];    
   
   
Bc= [  0   0    0   0   0   0   0   0   0   0
       1/m 0    0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   0
        0  1/m  0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   1/m 0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   1/m 0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   1/m 0   0   0   0   0  
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   1/m 0   0   0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   0   1/m 0   0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   1/m 0   0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   1/m 0
        0   0   0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0   0   1/m ];
% Bc= [  0   0   0    0   0   0   
%        1/m 0   0    0   0   0   
%         0   0  0    0   0   0   
%         0  1/m 0    0   0   0
%         0   0  0    0   0   0   
%         0   0  1/m  0   0   0
%         0   0   0   0   0   0
%         0   0   0   1/m 0   0
%         0   0   0   0   0   0   
%         0   0   0   0   1/m 0     
%         0   0   0   0   0   0   
%         0   0   0   0   0   1/m
%  ];


C=[ 1  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0  0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0  0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
    0  0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0
    0  0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0
    0  0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0
    0  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0
    0  0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0
    0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0
    0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0];
%   C=[ 1  0   0   0   0   0   0   0   0   0   0   0
%       0  0   1   0   0   0   0   0   0   0   0   0
%       0  0   0   0   1   0   0   0   0   0   0   0
%       0  0   0   0   0   0   1   0   0   0   0   0   
%       0  0   0   0   0   0   0   0   1   0   0   0   
%       0  0   0   0   0   0   0   0   0   0   1   0];


A=eye(20) + Ac*Ts;
B = Bc*Ts;
nu = size(B,2);
ny = size(C,1);
nx = size(B,1);

% 
% A_ext =  [A zeros(nx,ny);C*A eye(ny)];
% B_ext = [B;C*B]; 
% C_ext = [zeros(ny,nx),eye(ny)];


% Qi=10;%2600
% Ri=0.0002*Qi;%0.00017;
% Qi = 153;%100
% Ri = 0.000016; 
% Pi=100;
% Pi = 10;
Qi =3800;%100
Ri = 0.0001%0.00019;
Pi=500000;
N=50;

P = blkdiag(Pi,Pi,Pi,Pi,Pi,Pi,Pi,Pi,Pi,Pi);
Q = blkdiag(Qi,0.9*Qi,0.85*Qi,0.8*Qi,0.75*Qi,0.7*Qi, 0.65*Qi, 0.6*Qi, 0.55*Qi, 0.5*Qi);
R = blkdiag(Ri,0.9*Ri,0.85*Ri,0.8*Ri,0.75*Ri,0.7*Ri, 0.65*Ri, 0.6*Ri, 0.55*Ri, 0.5*Ri);

nk = 600;%250;
TU = 1:nk;
TX = 1:nk+1;
%TX = 0:T:(nk+1)*T - T;
Tref = 1:nk+N;
ref = 10 * square(0.0002*Tref, 0.79);
ref = [ref;ref;ref;ref;ref;ref;ref;ref;ref;ref];
%%
x10 = [1 0]';
x20 = [1 0]';
x30 = [1 0]';
x40 = [1 0]';
x50 = [1 0]';
x60 = [1 0]';
x70 = [1 0]';
x80 = [1 0]';
x90 = [1 0]';
x100 = [1 0]';
xd0 = [x10; x20; x30; x40;x50;x60;x70;x80;x90;x100];


[Fc,Gc,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(A,B,C,N,P,Q,R);

Fb = H*Fc;
Gb = H*Gc;
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^-1*St;
K = -Ky*Fb;

x0 = [xd0*0 ; C*xd0];
U = zeros(nu,nk);
U_d = zeros(nu,nk);
dU = zeros(nu,N);
dU_d = zeros(nu,N);
Xd(:,1) = xd0;
X(:,1) = x0;
Y(:,1) = C*xd0;
Xd(:,2) = xd0;
X(:,2) = x0;
Y(:,2) = C*xd0;


%compute constraints
u_max = fmax -(0.5*rho*area*Cd*ve^2);
u_min = fmin - (0.5*rho*area*Cd*ve^2);
U_max = kron(u_max,ones(N*nu,1));
U_min = kron(u_min,ones(N*nu,1));
M3 = kron(tril(ones(N)), eye(nu));
M4 = kron(ones(N,1), eye(nu));
Mu = [-M3;M3];

du_max = deltaFmax;
du_min = -deltaFmax;
DU_max = kron(du_max,ones(N*nu,1));
DU_min = kron(du_min,ones(N*nu,1));
Mdu = [-eye(N*nu);eye(N*nu)];
wdu = [-DU_min;DU_max];


M = [Mu;Mdu];
%M = Mu;
U(:,1) = 0;


for k = 2:nk
    
    Yb = reshape(ref(:,k:k+N),[],1);
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k);
    u_1 = U(:,k-1);
     wu = [-U_min + M4*u_1;U_max - M4*u_1];
     w = [wu;wdu];
    
    %w = wu;
    % centralized MPC
      [dUo,Jo,exitflag,output,lambda] = quadprog_v2(2*Rt,2*St*(Fb*xk-Yb),M,w);
    if exitflag<0
            
            error('Problems in the Optimization problem.');
            
            
    else
        flag = 1;
    end    
    
    %dUopt(:,:,k) = reshape( K*X(:,k)+Ky*Yb ,[],N);
    dUopt(:,:,k) = reshape( dUo ,nu,N); 
    duPlot(:,k) = dUopt(:,1,k) ;
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    %Xopt(:,:,k) = reshape( Fc*xk-Gc*(K*xk-Ky*Yb) ,6,N+1);
    U(:,k) = Uopt(:,1,k);
    
    Xd(:,k+1) = A*Xd(:,k) + B*U(:,k) + 0.005*rand;
    Y(:,k+1) = C*Xd(:,k+1);

    
    
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];
%MAX = u_max*ones(1,nk);

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
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstdarkblue);
plot(TX,Xd(3,:),'.--','Color',sstlightblue);
plot(TX,Xd(4,:),'-*','Color','magenta');
plot(TX,Xd(5,:),'.--','Color','red');
plot(TX,Xd(6,:),'-*','Color','yellow');
plot(TX,Xd(7,:),'.--','Color','black');
plot(TX,Xd(8,:),'-*','Color','cyan');
plot(TX,Xd(9,:),'.--','Color',sstgreen);
plot(TX,Xd(10,:),'-*','Color','green');
plot(TX,Xd(11,:),'.--','Color',sstgray);
plot(TX,Xd(12,:),'-*','Color',sstdarkgreen);
plot(TX,Xd(13,:),'.--','Color','#D95319');
plot(TX,Xd(14,:),'-*','Color','#EDB120');
plot(TX,Xd(15,:),'.--','Color','#7E2F8E');
plot(TX,Xd(16,:),'-*','Color','#77AC30');
plot(TX,Xd(17,:),'.--','Color','#4DBEEE');
plot(TX,Xd(18,:),'-*','Color','#A2142F');
plot(TX,Xd(19,:),'.--','Color','blue');
plot(TX,Xd(20,:),'-*','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r$$','car1 $$x_1$$','car1 $$x_2$$','car 2$$x_1$$','car 2 $$x_2$$','car 3$$x_1$$','car 3$$x_2$$','car 4$$x_1$$','car 4$$x_2$$','car 5$$x_1$$','car 5$$x_2$$','car 6$$x_1$$','car 6$$x_2$$' ,'car 7$$x_1$$','car 7$$x_2$$','car 8$$x_1$$','car 8$$x_2$$','car 9$$x_1$$','car 9$$x_2$$','car 10$$x_1$$','car 10$$x_2$$','Location','SouthEast');
title('State evolution');

figure(7103);
plot(TU,U(1,:),'s-','Color','magenta');
grid on;
hold on;
plot(TU,U(2,:),'d-','Color',sstdarkblue);
plot(TU,U(3,:),'d-','Color','red');
plot(TU,U(4,:),'d-','Color','cyan');
plot(TU,U(5,:),'d-','Color','black');
plot(TU,U(6,:),'d-','Color','green');
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ ','$$u_2$$ ', '$$u_3$$', '$$u_4$$', '$$u_5$$', '$$u_6$$');
title('Input');











   
   
   
   