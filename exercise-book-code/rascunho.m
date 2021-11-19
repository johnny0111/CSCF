clear all;

Ad = [0.3 0.6;0 1.2]; Bd = [0;1]; Cd = [1,0]; xd0 = [ -2;2];
P = 1;
Q = 1;
R = 4;
N = 3;
nxd = size(Bd ,1);
nu = size(Bd ,2);
ny = size(Cd ,1);
nx = nxd+ny;
[F,G,Qb ,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(Ad,Bd,Cd ,N,P,Q,R);
Fb = H*F; Gb = H*G;
Fdb = Hd*Fd; Gdb = Hd*Gd;
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^( -1)*St;
K = Rt^(-1)*St*Fb;
u_max = 0.8;
U_max = kron(u_max ,ones(N ,1));
M3 = tril(ones(N*nu));
M4 = ones(N*nu ,nu);
Mu = [-M3;M3];
M = Mu;
%simulate controlled system:
nk = 60; TU = 1:nk; TX = 1:nk+1; Tref = 1:nk+N; ref = 2*(Tref >=0);
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
    
Yb = ref(:,k:k+N)';
Dxdk = Xd(:,k)-Xd(:,k-1);
xk = [ Dxdk; Cd*Xd(:,k)];
u_1 = U(:,k-1);
wu = [U_max + M4*u_1;U_max - M4*u_1];
w = wu;
[dUo,Jo,exitflag,output,lambda] = quadprog(2*Rt ,2*St*(Fb*xk-Yb),M,w);
dUopt(:,:,k) = reshape( dUo ,nu,N);
Xopt(:,:,k) = reshape( F*xk-G*(K*xk -Ky*Yb) ,nx,N+1);
TUopt(k,:) = k:k+N-1; TXopt(k,:) = k:k+N;
Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
U(:,k) = Uopt(:,1,k);
Xd_opt(:,:,k) = reshape( Fd*Xd(:,k)+Gd*Uopt(:,:,k)' ,nxd ,N+1);
Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k);
X(:,k+1) = Xopt(:,2,k);

% get unconstrained response
Dxdk2 = Xd2(:,k)-Xd2(:,k-1);
xk2 = [ Dxdk2; Cd*Xd2(:,k)];
dUopt2 = reshape(-(K*xk2 -Ky*Yb) ,nu,N);
U2(:,k) = U2(:,k-1) + dUopt2 (: ,1);
Xd2(:,k+1) = Ad*Xd2(:,k) + Bd*U2(:,k);
end

figure (4401); plot(Xd(1,:),Xd(2,:),'bs-');
grid on; hold on; plot(Xd2(1,:),Xd2(2,:),'k+-');
for k = 1:nk, plot(Xd_opt(1,:,k),Xd_opt(2,:,k),'c.-'); end
plot(Xd(1,:),Xd(2,:),'bs-'); hold off; xlabel('x_1'); ylabel('x_2');
legend('const.','unconst.','pred.'); title('Phase plot');
figure (4402); plot(Tref ,ref ,'gs-'); grid on; hold on;
plot(TX ,Xd(1,:),'bs-'); plot(TX,Xd(2,:),'rd -');
plot(TX ,Xd2(1,:),'k+-'); plot(TX,Xd2(2,:),'mx-');
for k = 1:nk, plot(TXopt(k,:), Xd_opt(1,:,k),'c.-');
plot(TXopt(k,:), Xd_opt(2,:,k),'y.-'); end
plot(TX ,Xd(1,:),'bs-'); plot(TX,Xd(2,:),'rd -');
hold off; xlabel('k'); ylabel('x(k)');
legend('ref','x_1 con.','x_2 con.','x_1 unc.','x_2 unc.','x_1 pred.','x_2 pred.');
title('State evolution');
figure (4403); plot(TU ,U,'bs-'); grid on; hold on; plot(TU,U2 ,'k+-');
for k = 1:nk, plot(TUopt(k,:),Uopt(:,:,k),'c.-'); end
plot(TU ,U,'bs-'); hold off; xlabel('k'); ylabel('u(k)');
legend('const.','unc.','pred.'); title('Input ');