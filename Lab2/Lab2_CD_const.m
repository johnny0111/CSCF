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

m = 1500;
Cd = 0.3;
area = 1.8;
rho = 1.225;
g = 9.8065;
theta = 0.05;
fmax = 2800;
fmin = -2200;
deltaFmax = 200;
maxTheta = 0.1;
dDist = 10;
distMin = 3;
vRe = 20;
ve = 25;
we = 0;
T = 0.1;

N = 30;
Pi = 1;
Qi = 1000;
Ri = 0.001;

P = blkdiag(Pi,Pi);
Q = blkdiag(Qi,Qi);
R = blkdiag(Ri,Ri);

nk = 500;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = 10 * square(0.0002*Tref, 0.79);
ref = [ref; ref];
%%
%centralized
Ac = [0 -1 0 0; 0 -(rho*area*Cd*ve)/m 0 0; 0 1 0 -1; 0 0 0 -(rho*area*Cd*ve)/m];
Bc=[0 0; 1/m 0; 0 0; 0 1/m];
C=[1 0 0 0; 0 0 1 0];
A = eye(4) + Ac*T;
B = Bc*T;
x10 = [1 0]';
x20 = [1 0]';
xd0 = [x10; x20];


% compute centralized tracking controller
[Fc,Gc,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(A,B,C,N,P,Q,R);

Fb = H*Fc;
Gb = H*Gc;
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^-1*St;
K = -Ky*Fb;



nu = size(B,2);
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


%%compute constraints
u_max = fmax -(0.5*rho*area*Cd*ve^2);
u_min = -fmin - (0.5*rho*area*Cd*ve^2);
U_max = kron(u_max,ones(N*nu,1));
U_min = kron(u_min,ones(N*nu,1));
M3 = tril(ones(N*nu));
M4 = ones(N*nu,nu);
Mu = [-M3;M3];

du_max = deltaFmax;
du_min = -deltaFmax;
DU_max = kron(du_max,ones(N*nu,1));
DU_min = kron(du_min,ones(N*nu,1));
Mdu = [-eye(N*nu);eye(N*nu)];
wdu = [-DU_min;DU_max];

pr_min = -7 ;
pr_max = 93;
Y_min = kron(pr_min,ones(nu*(N+1),1));
Y_max = kron(pr_max, ones(nu*(N+1),1));
My = [-Gb; Gb];

M = [Mu;Mdu;My];

U(:,1) = 0;
%%


for k = 2:nk
    
    % compute initial condition and current reference sequence
%     Yb = [ref(:,k:k+N)';ref(:,k:k+N)'];
%     Yb1 = ref(1,k:k+N)';
%     Yb2 = ref(2,k:k+N)';
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k);
    u_1 = U(:,k-1);
    wu = [U_max + M4*u_1;U_max - M4*u_1];
    wy = [-Y_min + Fb*xk; Y_max - Fb*xk];
    w = [wu;wdu;wy];
    % centralized MPC
      [dUo,Jo,exitflag,output,lambda] = quadprog_v2(2*Rt,2*St*(Fb*xk-Yb),M,w);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end    
    
    %dUopt(:,:,k) = reshape( K*X(:,k)+Ky*Yb ,[],N);
    dUopt(:,:,k) = reshape( dUo ,nu,N); 
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    Xopt(:,:,k) = reshape( Fc*xk-Gc*(K*xk-Ky*Yb) ,6,N+1);
    uk = Uopt(:,1,k);
    U(:,k) = uk;
    
    Xd(:,k+1) = A*Xd(:,k) + B*U(:,k) ;
    Y(:,k+1) = C*Xd(:,k+1);

    
    % Decentralized MPC
%     Dxdk_d = Xd_d(:,k) - Xd_d(:,k-1);
%     X_d(:,k) = [Dxdk_d; C1*Xd_d(1:2,k); C2*Xd_d(3:4,k)];
    
%     dU1opt(:,:,k) = reshape( K1*[X_d(1:2,k); X_d(5,k)]+Ky1*Yb1 ,[],N);
%     %X1opt(:,:,k) = reshape( F1*[Xd_d(1:2,k); X_d(5,k)] +G1*(K1*[Xd_d(1:2,k); X_d(5,k)]+Ky1*Yb1) ,nx1,N+1);
%     dU2opt(:,:,k) = reshape( K2*[X_d(3:4,k); X_d(6,k)]+Ky2*Yb2 ,[],N);
%     %X2opt(:,:,k) = reshape( F2*Xd_d(2,k)+G2*(K2*Xd_d(2,k)+Ky2*Yb2) ,nx2,N+1);
%     dUd(:,k) = [dU1opt(:,1,k);dU2opt(:,1,k)];   
%     
%     Uopt_d(:,:,k) = U_d(:,k-1) + dUd(:,k);
%     % joint system simulation for decentralized MPC 
%     Xd_d(:,k+1) = A_d*Xd_d(:,k) + B_d*Uopt_d(:,k);
%     % compute auxiliary variables for visualization:
    
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];

%%
figure(7101);
%plot(ref(1,:),ref(2,:),'k+-');
grid on;
hold on;
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
plot(Xd(3,:), Xd(4,:), 's-','Color','magenta');
%plot(Xd_d(1,:),Xd_d(2,:),'o-','Color',sstgreen);
%plot(Xd_d(3,:),Xd_d(4,:),'o-','Color','red');

% plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
% plot(Xd_d(1,:),Xd_d(2,:),'o-','Color',sstgreen);
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('car 1 cent.','car 2 cent','car 1 dec.','car2 dec.');
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
% plot(TX,Xd_d(1,:),'o-','Color',sstgreen);
% plot(TX,Xd_d(2,:),'+-','Color',sstdarkgreen);
% plot(TX,Xd_d(3,:),'o-','Color','black');
% plot(TX,Xd_d(4,:),'+-','Color','yellow');
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','car1 $$x_1$$ cen.','car1 $$x_2$$ cen.','car 2$$x_1$$ cen.','car 2 $$x_2$$ cen.','car 1$$x_1$$ decen.','car 1$$x_2$$  deccen.','car 2$$x_1$$ decen.','car 2$$x_2$$  deccen.','Location','SouthEast');
title('State evolution');

figure(7103);
plot(TU,U(1,:),'s-','Color',sstblue);
grid on;
hold on;
plot(TU,U(2,:),'d-','Color',sstdarkblue);
% plot(TU,dUd(1,:),'o-','Color',sstgreen);
% plot(TU,dUd(2,:),'+-','Color',sstdarkgreen);
plot(TU,U(1,:),'s-','Color',sstblue);
plot(TU,U(2,:),'d-','Color',sstdarkblue);
% plot(TU,dUd(1,:),'o-','Color',sstgreen);
% plot(TU,dUd(2,:),'+-','Color',sstdarkgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.','$$u_1$$ decent.','$$u_2$$ decent.');
title('Input');









