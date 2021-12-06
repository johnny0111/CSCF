%%
function flag = Lab2_CD_const(Ri,Qi)
%%
% function flag = Lab2_CD_const(Ri,Qi)
%clear all;
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
vRe = 25;
ve = 25;
we = 0;
T = 0.1;

N=45;

% Qi=40;%2600
% Ri=0.0002*Qi;%0.00017;


Pi=100;


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
%centralized
Ac = [0 -1 0 0; 0 -(rho*area*Cd*ve)/m 0 0; 0 1 0 -1; 0 0 0 -(rho*area*Cd*ve)/m];
Bc=[0 0; 1/m 0; 0 0; 0 1/m];
C=[1 0 0 0; 0 0 1 0];
A = eye(4) + Ac*T;
B = Bc*T;
x10 = [5 2]';
x20 = [5 2]';
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

pr_min = -7 ;
pr_max = 93;
Y_min = kron(pr_min,ones(nu*(N+1),1));
Y_max = kron(pr_max, ones(nu*(N+1),1));
My = [-Gb; Gb];

M = [Mu;Mdu];
%M = Mu;
U(:,1) = 0;
%%


for k = 2:nk
    
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k);
    u_1 = U(:,k-1);
    wu = [-U_min + M4*u_1;U_max - M4*u_1];
    wy = [-Y_min + Fb*xk; Y_max - Fb*xk];
    w = [wu;wdu];
    
    %w = wu;
    % centralized MPC
      [dUo,Jo,exitflag,output,lambda] = quadprog_v2(2*Rt,2*St*(Fb*xk-Yb),M,w);
    if exitflag~=1
            flag = -1
            return;
    else
        flag = 1;
    end    
    
    %dUopt(:,:,k) = reshape( K*X(:,k)+Ky*Yb ,[],N);
    dUopt(:,:,k) = reshape( dUo ,nu,N); 
    duPlot(:,k) = dUopt(:,1,k) ;
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    Xopt(:,:,k) = reshape( Fc*xk-Gc*(K*xk-Ky*Yb) ,6,N+1);
    U(:,k) = Uopt(:,1,k);
    
    Xd(:,k+1) = A*Xd(:,k) + B*U(:,k) ;
    Y(:,k+1) = C*Xd(:,k+1);

    
    
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];
MAX = u_max*ones(1,nk);

%%


figure(7101);

grid on;
hold on;
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
plot(Xd(3,:), Xd(4,:), 's-','Color','magenta');

hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('car1 cent.','car2 cent','car1 dec.','car2 dec.');
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
legend('$$r_1$$','$$r_2$$','car1 $$x_1$$ cen.','car1 $$x_2$$ cen.','car 2$$x_1$$ cen.','car 2 $$x_2$$ cen.','car 1$$x_1$$ decen.','car 1$$x_2$$  deccen.','car 2$$x_1$$ decen.','car 2$$x_2$$  deccen.','Location','SouthEast');
title('State evolution');

figure(7103);
plot(TU,U(1,:),'s-','Color','magenta');
grid on;
hold on;
plot(TU,U(2,:),'d-','Color',sstdarkblue);
plot(TU, duPlot(1,:), 'k+--','Color', sstdarkgreen);
plot(TU, duPlot(2,:), '-*','Color', sstlightgreen);
plot(TU, MAX);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.', 'du1', 'du2');
title('Input');


%end









end






