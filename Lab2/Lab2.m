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

%%
Ac = [0 -1 0 0; 0 -(rho*area*Cd*ve)/m 0 0; 0 1 0 -1; 0 0 0 -(rho*area*Cd*ve)/m];
Bc=[0 0; 1/m 0; 0 0; 0 1/m];
C=[1 0 0 0; 0 0 1 0];
A = eye(4) + Ac*T;
B = Bc*T;
x10 = [1 0]';
x20 = [1 0]';
xd0 = [x10; x20];

N = 20;
Pi = 1;
Qi = 1000;
Ri = 0.001;


% compute centralized tracking controller
P = blkdiag(Pi,Pi);
Q = blkdiag(Qi,Qi);
R = blkdiag(Ri,Ri);
[Fc,Gc,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(A,B,C,N,P,Q,R);

Fb = H*Fc;
Gb = H*Gc;
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^-1*St;
K = -Ky*Fb;


nk = 500;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
%ref(1:length(Tref)) = -10*ones(); % step reference of amplitude 1
%ref(351:length(Tref)) = -10*ones();
ref = 10 * square(0.0002*Tref, 0.79);
%ref = 10 * square(0.0033*Tref,7.5);
%dist_x1 = 0*0.5*ones(size(ref)).*(Tref>=20);
nu = size(B,2);
 x0 = [xd0*0 ; C*xd0];
 U = zeros(nu,nk);
dU = zeros(nu,N);
Xd(:,1) = xd0;
X(:,1) = x0;
Y(:,1) = C*xd0;
Xd(:,2) = xd0;
X(:,2) = x0;
Y(:,2) = C*xd0;



for k = 2:nk
    
    % compute initial condition and current reference sequence
    Yb = [ref(:,k:k+N)';ref(:,k:k+N)'];
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; C*Xd(:,k)];
    xk = X(:,k);
    
    % centralized MPC
    dUopt(:,:,k) = reshape( K*X(:,k)+Ky*Yb ,[],N);
    
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    Xopt(:,:,k) = reshape( Fc*xk-Gc*(K*xk-Ky*Yb) ,6,N+1);
    uk = Uopt(:,1,k);
    U(:,k) = uk;
    
    Xd(:,k+1) = A*Xd(:,k) + B*U(:,k) ;
    Y(:,k+1) = C*Xd(:,k+1);

    % compute auxiliary variables for visualization:


    
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; C*Xd(:,k+1)];

figure(7101);
plot(ref(1,:),ref(1,:),'k+-');
grid on;
hold on;
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
%plot(Xd(1,:),Xd(2,:),'o-','Color',sstgreen);
for k = 1:nk
    plot(Xopt(1,:,k),Xopt(2,:,k),'.--','Color',sstlightblue);
    %plot(X1opt(1,:,k),X2opt(1,:,k),'.-.','Color',sstlightgray);
end
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
%plot(Xd(1,:),Xd(2,:),'o-','Color',sstgreen);
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('ref.','cent.','decent.','pred. cent.','pred. dec.');
title('Phase plot');

figure(7102);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(Tref,ref(1,:),'k+--');
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstdarkblue);
%plot(TX,Xd(1,:),'o-','Color',sstgreen);
%plot(TX,Xd(2,:),'+-','Color',sstdarkgreen);
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-.','Color',sstlightgray);
end
plot(Tref,ref(1,:),'k+-');
plot(Tref,ref(1,:),'k+--');
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstdarkblue);
%plot(TX,Xd(1,:),'o-','Color',sstgreen);
%plot(TX,Xd(2,:),'+-','Color',sstdarkgreen);
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
%plot(TU,Ud(1,:),'o-','Color',sstgreen);
%plot(TU,Ud(2,:),'+-','Color',sstdarkgreen);
for k = 1:nk
    plot(TUopt(k,:),Uopt(1,:,k),'.--','Color',sstlightblue);
    plot(TUopt(k,:),Uopt(2,:,k),'.-.','Color',sstlightgray);
end
plot(TU,U(1,:),'s-','Color',sstblue);
plot(TU,U(2,:),'d-','Color',sstdarkblue);
%plot(TU,Ud(1,:),'o-','Color',sstgreen);
%plot(TU,Ud(2,:),'+-','Color',sstdarkgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.','$$u_1$$ decent.','$$u_2$$ decent.');
title('Input');



