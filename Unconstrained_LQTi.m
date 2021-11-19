clear all, clc

sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;

%% Unconstrained LQTi

clear X U Xopt Uopt TX TU TUopt TXopt;

Ad = [1 1;0 1];
Bd = [0;1];
Cd = [1 0];
N = 3;
P = 1;
Q = 1;
R = 4;
Rb = R*eye(N);
xd0 = [-4.5;2];
nxd = size(Bd,1);
nu = size(Bd,2);
ny = size(Cd,1);
nx = nxd+ny;

% compute batch matrices
[F,G,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
Fb = H*F;
Gb = H*G;
Fdb = Hd*Fd;
Gdb = Hd*Gd;
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^(-1)*St,
K = Rt^(-1)*St*Fb,

%simulate controlled system:
nk = 60;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = (Tref>=30); % step reference of amplitude 1, starting at k=30
x0 = [xd0*0 ; Cd*xd0];
U = zeros(nu,nk);
U2 = zeros(nu,nk);
Xd(:,1) = xd0;
Xd2(:,1) = xd0;
X(:,1) = x0;
Xd(:,2) = xd0;
Xd2(:,2) = xd0;
X(:,2) = x0;
for k = 2:nk    
    % compute initial conditions
    Yb = ref(:,k:k+N)';    
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; Cd*Xd(:,k)];
    xk = X(:,k);
    u_1 = U(:,k-1);
    Dxdk2 = Xd2(:,k)-Xd2(:,k-1);
    xk2 = [ Dxdk2; Cd*Xd2(:,k)];
    
    % compute optimal control sequence:
    dUopt(:,:,k) = reshape( -K*X(:,k) + Ky*Yb ,nu,N);
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);

    % set MPC control policy:
    U(:,k) = Uopt(:,1,k);
    
    % simulate system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k);
    
    % compute unconstrained optimal sequence and MPC policy
    dUopt2(:,:,k) = reshape(-(K*xk2-Ky*Yb) ,nu,N);
    Uopt2(:,:,k) = U2(:,k-1) + dUopt2(:,:,k);
    U2(:,k) = Uopt2(:,1,k);
    
    % simulate system 2 (for comparison):
    Xd2(:,k+1) = Ad*Xd2(:,k) + Bd*U2(:,k);
    
    % compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( F*xk-G*(K*xk-Ky*Yb) ,nx,N+1);
    Xd_opt(:,:,k) = reshape( Fd*Xd(:,k)+Gd*Uopt(:,:,k)' ,nxd,N+1);
    Xd_opt2(:,:,k) = reshape( Fd*Xd2(:,k)+Gd*Uopt2(:,:,k)' ,nxd,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; Cd*Xd(:,k+1)];
    
figure(4501);
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
grid on;
hold on;
plot(Xd2(1,:),Xd2(2,:),'o-','Color',sstgray);
for k = 1:nk
    plot(Xd_opt(1,:,k),Xd_opt(2,:,k),'.--','Color',sstlightblue);
    plot(Xd_opt2(1,:,k),Xd_opt2(2,:,k),'.-.','Color',sstlightgray);
end
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
plot(Xd2(1,:),Xd2(2,:),'o-','Color',sstgray);
hold off;
xlabel('x_1');
ylabel('x_2');
legend('const.','unconst.','const. pred.','unconst. pred.');
title('Phase plot');

figure(4502);
plot(Tref,ref,'kx-');
grid on;
hold on;
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstgreen);
plot(TX,Xd2(1,:),'o-','Color',sstgray);
plot(TX,Xd2(2,:),'x-','Color',sstlightgray);
for k = 1:nk
    plot(TXopt(k,:),Xd_opt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xd_opt(2,:,k),'.-.','Color',sstlightgreen);
end
plot(TX,Xd(1,:),'s-','Color',sstblue);
plot(TX,Xd(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
legend('ref','x_1 const.','x_2 const.','x_1 unc.','x_2 unc.','x_1 pred.','x_2 pred.','Location','SouthEast');
title('State evolution');

figure(4503);
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
plot(TU,U2,'o-','Color',sstgray);
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.--','Color',sstlightblue);
    plot(TUopt(k,:),Uopt2(:,:,k),'.-.','Color',sstlightgray);
end
plot(TU,U,'s-','Color',sstblue);
hold off;
xlabel('k');
ylabel('u(k)');
legend('const.','unc.','pred.');
title('Input');