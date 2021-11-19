%%
clear all;
close all;
sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;
clc;
%Problem parameters
m = 1500;
Cd = 0.3;
area = 1.8;
rho = 1.225;
g = 9.8065;
theta = 0.05;
fmax = 2600;
fmin = -2000;
deltaFmax = 50;
maxTheta = 0.1;
dDist = 10;
distMin = 5;
vRe = 20;
ve = 20;
we = 0;
T = 0.1;
%N=50;
N=70;
%R=0.1;
R=0.0000001
Q = 100;
%Q=0.5;
P = 100;
%P=300;
%%
%%%%%%%%%%%%%%%%%%% System modulation and descretization
Ac = [0 -1; 0 (rho*area*Cd*ve) / m];
Bc=[0; 1/m];
Bdist = [1 0 0;0 -g -(rho*area*Cd*2*ve)/(2*m)];
C=[1 0];
D = 0;
Ad = eye(2) + Ac*T;
Bd = Bc*T;
Bdistd = Bdist * T;
lambda = eig(Ad);

%simulation
sys = ss(Ad,Bd,eye(2),D,T);
xd0 = [-0.2;0.3];
nxd = size(Bd,1);
nu = size(Bd,2);
ny = size(C,1);
t = 0:T:500;
u = 30*ones(length(t),1);
%lsim(sys,u,t,xd0)
grid on
%%
%%%%%%%%%%%%%%%%%%%% Extended model analysis and unc. LQTi
AExt = [Ad zeros(2,1) ; C*Ad 1];
BExt = [Bd; C*Bd]
CExt = [zeros(1,2) 1];
nx = size(BExt,1);
lambdaAug = eig(AExt);

[F,G,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(Ad,Bd,C,N,P,Q,R);
Fb = H*F;
Gb = H*G;
Fdb = Hd*Fd;
Gdb = Hd*Gd;

% compute final cost function matrices and control gains
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^(-1)*St;
K  = Rt^(-1)*St*Fb;

%simulate controlled system:
nk = 500;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
%ref(1:length(Tref)) = -10*ones(); % step reference of amplitude 1
%ref(351:length(Tref)) = -10*ones();
ref = 10 * square(0.0002*Tref, 0.79);
%ref = 10 * square(0.0033*Tref,7.5);
%dist_x1 = 0*0.5*ones(size(ref)).*(Tref>=20);
x0 = [xd0*0 ; C*xd0];
U = zeros(nu,nk);
dU = zeros(nu,N);
Xd(:,1) = xd0;
X(:,1) = x0;
Y(:,1) = C*xd0;
Xd(:,2) = xd0;
X(:,2) = x0;
Y(:,2) = C*xd0;
%%
for k = 2:nk
    
    % compute initial condition and current reference sequence
    Yb = ref(:,k:k+N)';
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k);
    
    % compute optimal incremental control sequence:
    dUopt(:,:,k) = reshape(-(K*xk-Ky*Yb) ,nu,N);

    % set MPC control policy and simulate system:
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    uk = Uopt(:,1,k);
    U(:,k) = uk;
%     if U(:,k) > (2600 -(0.5*rho*area*Cd*ve^2))
%         U(:,k) = 2600 -(0.5*rho*area*Cd*ve^2);
%     end
%     if U(:,k) - U(:,k-1) > 50
%        U(:,k) = U(:,k-1) + 50; 
%     end
    % simulate original system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k) ;
    Y(:,k+1) = C*Xd(:,k+1);
    if(Xd(1,k+1) < -5)
       Xd(1,k+1) = -5; 
    end
    % compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( F*xk-G*(K*xk-Ky*Yb) ,nx,N+1);
    Yopt(:,:,k) = reshape( Fb*xk-Gb*(K*xk-Ky*Yb) ,ny,N+1);
    Xd_opt(:,:,k) = reshape( Fd*Xd(:,k)+Gd*Uopt(:,:,k)' ,nxd,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];
    %%
figure(3401);
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
grid on;
hold on;
for k = 1:nk
    plot(Xd_opt(1,:,k),Xd_opt(2,:,k),'.-','Color',sstlightblue);
end
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
hold off;
xlabel('x_1');
ylabel('x_2');
legend('trajectory','predictions');
title('Phase plot');

figure(3402);
plot(Tref,ref,'gs-');
grid on;
hold on;
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstgreen);
for k = 1:nk
    plot(TXopt(k,:),Xd_opt(1,:,k),'.-','Color',sstlightblue);
    plot(TXopt(k,:),Xd_opt(2,:,k),'.-','Color',sstlightgreen);
end
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
legend('ref','x_1','x_2','x_1 measured','x_1 predictions','x_2 predictions','Location','SouthEast');
title('State evolution');

figure(3403);
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.-','Color',sstlightblue);
end
plot(TU,U,'s-','Color',sstblue);
hold off;
xlabel('k');
ylabel('u(k)');
legend('input','input predictions');
title('Input');

figure(3404);
plot(TX,X(1,:),'s-','Color',sstblue);
grid on;
hold on;
plot(TX,X(2,:),'d-','Color',sstgreen);
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.-','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-','Color',sstlightgreen);
end
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
legend('\Delta x_1','\Delta x_2','\Delta x_1 predictions','\Delta x_2 predictions','Location','SouthEast');
title('\Delta State evolution');

%%

% compute constraints matrices:
u_max = 2600 -(0.5*rho*area*Cd*ve^2);
u_min = -2000 - (0.5*rho*area*Cd*ve^2);
du_max = 50;
du_min = -50;
pr_min = -5 ;
pr_max = 90 ;
U2 = zeros(nu,nk);

U_max = kron(u_max,ones(N,1));
U_min = kron(u_min,ones(N,1));
DU_max = kron(du_max,ones(N,1));
DU_min = kron(du_min,ones(N,1));
Y_min = kron(pr_min,ones(N+1,1));
Y_max = kron(pr_max, ones(N+1,1));

M3 = tril(ones(N*nu));
M4 = ones(N*nu,nu);
Mu = [-M3;M3];
Mdu = [-eye(N);eye(N)];
My = [-Gb; Gb];
M = [Mu;My;Mdu];

wdu = [-DU_min;DU_max];

for k = 2:nk
    
    % compute initial conditions
    Yb = ref(:,k:k+N)';    
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k);
    u_1 = U(:,k-1);
    wu = [U_max + M4*u_1;U_max - M4*u_1];
    wy = [-Y_min + Fb*xk; Y_max - Fb*xk];
    w = [wu;wy;wdu];
    %Dxdk2 = Xd2(:,k)-Xd2(:,k-1);
    %xk2 = [ Dxdk2; Cd*Xd2(:,k)];
    
    % compute constrained optimal incremental control sequence and MPC policy
    [dUo,Jo,exitflag,output,lambda] = quadprog_v2(2*Rt,2*St*(Fb*xk-Yb),M,w);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end    
    dUopt(:,:,k) = reshape( dUo ,nu,N);    
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    U(:,k) = Uopt(:,1,k);
    
    % simulate original system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k);
    
%     % compute unconstrained optimal sequence and MPC policy
%     dUopt2(:,:,k) = reshape(-(K*xk2-Ky*Yb) ,nu,N);
%     Uopt2(:,:,k) = U2(:,k-1) + dUopt2(:,:,k);
%     U2(:,k) = Uopt2(:,1,k);
%     
%     % simulate system 2 (for comparison):
%     Xd2(:,k+1) = Ad*Xd2(:,k) + Bd*U2(:,k);
    
    % compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( F*xk-G*(K*xk-Ky*Yb) ,nx,N+1);
    Xd_opt(:,:,k) = reshape( Fd*Xd(:,k)+Gd*Uopt(:,:,k)' ,nxd,N+1);
    %Xd_opt2(:,:,k) = reshape( Fd*Xd2(:,k)+Gd*Uopt2(:,:,k)' ,nxd,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];
    
figure(4401);
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
grid on;
hold on;
%plot(Xd2(1,:),Xd2(2,:),'o-','Color',sstgray);
for k = 1:nk
    plot(Xd_opt(1,:,k),Xd_opt(2,:,k),'.--','Color',sstlightblue);
    %plot(Xd_opt2(1,:,k),Xd_opt2(2,:,k),'.-.','Color',sstlightgray);
end
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
%plot(Xd2(1,:),Xd2(2,:),'o-','Color',sstgray);
hold off;
xlabel('x_1');
ylabel('x_2');
%legend('const.','unconst.','const. pred.','unconst. pred.');
legend('const.','const. pred.');
title('Phase plot');

figure(4402);
plot(Tref,ref,'kx-');
grid on;
hold on;
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstgreen);
%plot(TX,Xd2(1,:),'o-','Color',sstgray);
%plot(TX,Xd2(2,:),'x-','Color',sstlightgray);
for k = 1:nk
    plot(TXopt(k,:),Xd_opt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xd_opt(2,:,k),'.-.','Color',sstlightgreen);
end
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
legend('ref','x_1 const.','x_2 const.','x_1 pred.','x_2 pred.','Location','SouthEast');
title('State evolution');

figure(4403);
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
plot(TU,U2,'o-','Color',sstgray);
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.--','Color',sstlightblue);
    %plot(TUopt(k,:),Uopt2(:,:,k),'.-.','Color',sstlightgray);
end
plot(TU,U,'s-','Color',sstblue);
hold off;
xlabel('k');
ylabel('u(k)');
legend('const.','pred.');
title('Input');






