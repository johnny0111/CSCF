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
T = 0.1;

N=75;

% Qi=10;%2600
% Ri=0.0002*Qi;%0.00017;
Qi = 153;%100
Ri = 0.000016;
Pi=100;
% Pi = 10;
Qi =3800;%100
Ri = 0.0000001%0.00019;
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
A1c = [0 -1; 0 -(rho*area*Cd*ve)/m];
A2c = [0 -1; 0 -(rho*area*Cd*ve)/m];
B1c = [0; 1/m];
B2c = [0; 1/m];
C1 = [1 0];
C2 = [1 0];
A1 = eye(2) + A1c*T;
A2 = eye(2) + A2c*T;
B1 = B1c*T;
B2 = B2c*T;
nx1 = size(B1,1);
nx2 = size(B2,1);
nu1 = size(B1,2);
nu2 = size(B2,2);
A_d = blkdiag(A1,A2);
B_d = blkdiag(B1,B2);
C_d = blkdiag(C1,C2);

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


x0_d = [xd0*0 ; C1*x10;C2*x20];
Xd_d(:,1) = xd0;
X_d(:,1) = x0_d;
Xd_d(:,2) = xd0;
X_d(:,2) = x0_d;


%% compute constraints
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

% Decentralized states
x_d0 = [xd0*0 ; C*xd0];
Xd_d(:,1) = xd0;
X_d(:,1)=x_d0;
Xd_d(:,2) = xd0;
X_d(:,2)=x_d0;

% VER LAB2_V2(?) e 2.7!!

for k = 2:nk
    
    %?
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
    
    %?
    Dxdk1 = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k); 
    
    %?
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k);
    
    u1_1 = U1(:,k-1);
    u2_1 = U2(:,k-1);
    
    % Decentralized MPC
      [dUo1,Jo1,exitflag,output1,lambda1] = quadprog_v2(2*Rt1,2*St1*(Fb1*x1k-Yb1),M,w);
    if exitflag~=1
            k
            error('Problems in the Optimization problem (1)');          
    else
        flag = 1;
    end    
    
          [dUo2,Jo2,exitflag,output2,lambda2] = quadprog_v2(2*Rt2,2*St2*(Fb2*x2k-Yb2),M,w);
    if exitflag~=1
            k
            error('Problems in the Optimization problem (2)');         
    else
        flag = 1;
    end 
      
    % 1
    dU1opt(:,:,k) = reshape( dUo1 ,nu,N);
    du1Plot(:,k) = dU1opt(:,1,k) ;
    U1opt(:,:,k) = U1(:,k-1) + dU1opt(:,:,k);
    X1opt(:,:,k) = reshape( Fc1*xk1-Gc1*(K1*xk1-Ky1*Yb1) ,6,N+1);
    U1(:,k) = U1opt(:,1,k);
    
    % 2
    dU2opt(:,:,k) = reshape( dUo2 ,nu,N);
    du2Plot(:,k) = dU2opt(:,1,k) ;
    U2opt(:,:,k) = U2(:,k-1) + dU2opt(:,:,k);
    X2opt(:,:,k) = reshape( Fc2*xk2-Gc2*(K2*xk2-Ky2*Yb2) ,6,N+1);
    U2(:,k) = U2opt(:,1,k);
    
    % Joint system
    Ud(:,k) = [U1opt(:,1,k);U2opt(:,1,k)];   
    Xd(:,k+1) = A*Xd(:,k) + B*Ud(:,k) + 0.005*rand;
    
    % Visualization
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];
MAX = u_max*ones(1,nk);

%%


figure(7301);
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

figure(7302);
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

figure(7303);
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


%end