%% Exercise book solutions (CPCS)
%  Module 2: Linear model prdictive control
%
% Bruno Guerreiro (bj.guerreiro@fct.unl.pt)

clear all;
sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;

%% Exercise 2.1: Double integrator MPC regulator 

clear X U Xopt Uopt TX TU TUopt TXopt;

A = [1 1;0 1];
B = [0;1];
N = 3;
P = eye(2);
Q = eye(2);
R = 10;
x0 = [-4.5;2];
nx = size(B,1);
nu = size(B,2);

%compute final cost function matrices and control gain
[F,G,Qb,Rb] = GetBatchXMatrices(A,B,[],N,P,Q,R);
Qt = F'*Qb*F,
Rt = G'*Qb*G + Rb,
St = G'*Qb*F,
K = Rt^-1*St;

%simulate controlled system:
X(:,1) = x0;
nk = 30;
TU = 1:nk;
TX = 1:nk+1;
for k = 1:nk
    
    % compute optimal control sequence:
    
    Uopt(:,:,k) = reshape( -K*X(:,k) ,nu,N);

    % set MPC control policy:
    U(:,k) = Uopt(:,1,k);
    
    % simulate system:
    X(:,k+1) = A*X(:,k) + B*U(:,k);
    
    % compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( (F-G*K)*X(:,k) ,nx,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
end

% visualization:
figure(3101);
plot(X(1,:),X(2,:),'bs-');
grid on;
hold on;
for k = 1:nk
    plot(Xopt(1,:,k),Xopt(2,:,k),'.-','Color',sstlightblue);
end
plot(X(1,:),X(2,:),'s-','Color',sstblue);
hold off;
xlabel('x_1');
ylabel('x_2');
legend('trajectory','predictions');
title('Phase plot');

figure(3102);
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
grid on;
hold on;
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.-','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-','Color',sstlightgreen);
end
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
legend('x_1','x_2','x_1 predictions','x_2 predictions','Location','SouthEast');
title('State evolution');

figure(3103);
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.-','Color',sstlightblue);
end
plot(TU,U,'bs-');
hold off;
xlabel('k');
ylabel('u(k)');
legend('input','input predictions');
title('Input');

%% Exercise 2.2: no solution provided.

%% Exercise 2.3: no solution provided.  LQT unconstraint
clear X U Xopt Uopt TX TU TUopt TXopt;

A = [0.4 1;0 3];
B = [0;2];
C = [1 0];
N = 3;
Yb = ones(N+1,1);
P = 1;
Q = 1;
R = 2;
x0 = [2;2];
nx = size(B,1);
nu = size(B,2);

%compute final cost function matrices and control gain
[F,G,Qb,Rb,H] = GetBatchXMatrices(A,B,C,N,P,Q,R);
%Qt = F'*Qb*F,
Fb = H*F,
Gb = H*G,
Rt = Gb'*Qb*Gb + Rb,
St = Gb'*Qb,
Ky = Rt^-1*St,
K = Ky*Fb

X(:,1) = x0;
nk = 30;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = (Tref>=1);
for k = 1:nk
    
    xk = X(:,k);
    % compute optimal control sequence:
    Yb = ref(:,k:k+N)';
    Uopt(:,:,k) = reshape(-K*xk+Ky*Yb ,nu,N);

    % set MPC control policy:
    uk = Uopt(:,1,k);
    U(:,k) = uk;
    % simulate system:
    X(:,k+1) = A*xk + B*U(:,k);

    % compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( (F-G*K)*xk + G*Ky*Yb,nx,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
end

% visualization:
figure();
plot(X(1,:),X(2,:),'bs-');
grid on;
hold on;
for k = 1:nk
    plot(Xopt(1,:,k),Xopt(2,:,k),'.-','Color',sstlightblue);
end
plot(X(1,:),X(2,:),'s-','Color',sstblue);
hold off;
xlabel('x_1');
ylabel('x_2');
legend('trajectory','predictions');
title('Phase plot');

figure();
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
grid on;
hold on;
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.-','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-','Color',sstlightgreen);
end
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
legend('x_1','x_2','x_1 predictions','x_2 predictions','Location','SouthEast');
title('State evolution');

figure();
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.-','Color',sstlightblue);
end
plot(TU,U,'bs-');
hold off;
xlabel('k');
ylabel('u(k)');
legend('input','input predictions');
title('Input');
%% Exercise 2.4: Double integrator tracking MPC with measurement errors

clear X U Xopt Uopt TX TU TXopt TUopt Xd ;

Ad = [1 1;0 1];
Bd = [0;1];
Cd = [1,0];
P = 1; %notice P has dimensions compatible with the output
Q = 1; %notice Q has dimensions compatible with the output
R = 4;
N = 3;
xd0 = [-4.5;2];
nxd = size(Bd,1);
nu = size(Bd,2);
ny = size(Cd,1);

A = [Ad , zeros(nxd,ny) ; Cd*Ad , eye(ny)];
B = [Bd ; Cd*Bd];
C = [zeros(ny,nxd),eye(ny)];
nx = size(B,1);

% compute matrices
[F,G,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
Fb = H*F;
Gb = H*G;
Fdb = Hd*Fd;
Gdb = Hd*Gd;

% compute final cost function matrices and control gains
Rt = Gb'*Qb*Gb + Rb,
St = Gb'*Qb,
Ky = Rt^(-1)*St,
K  = Rt^(-1)*St*Fb,

%simulate controlled system:
nk = 60;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = (Tref>=40); % step reference of amplitude 1
dist_x1 = 0*0.5*ones(size(ref)).*(Tref>=20);
x0 = [xd0*0 ; Cd*xd0];
U = zeros(nu,nk);
dU = zeros(nu,N);
Xd(:,1) = xd0;
X(:,1) = x0;
Y(:,1) = Cd*xd0;
Xd(:,2) = xd0;
X(:,2) = x0;
Y(:,2) = Cd*xd0;
for k = 2:nk
    
    % compute initial condition and current reference sequence
    Yb = ref(:,k:k+N)';
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; Cd*Xd(:,k)];
    xk = X(:,k);
    
    % compute optimal incremental control sequence:
    dUopt(:,:,k) = reshape(-(K*xk-Ky*Yb) ,nu,N);

    % set MPC control policy and simulate system:
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    uk = Uopt(:,1,k);
    U(:,k) = uk;
    
    % simulate original system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k) + [dist_x1(:,k);0];
    Y(:,k+1) = Cd*Xd(:,k+1);
    
    % compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( F*xk-G*(K*xk-Ky*Yb) ,nx,N+1);
    Yopt(:,:,k) = reshape( Fb*xk-Gb*(K*xk-Ky*Yb) ,ny,N+1);
    Xd_opt(:,:,k) = reshape( Fd*Xd(:,k)+Gd*Uopt(:,:,k)' ,nxd,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; Cd*Xd(:,k+1)];
    
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

%% Exercise 2.5: Constrained minimization

Q = eye(3);
s = [-1;-2;1];
M = [1 -3 -2; 2 1 2; 1 -2 2];
w = [-1;1;2];

% check for inactive constraints:
det(M),
lbd = -( M*2*Q^(-1)*M' )^(-1)*(4*M*Q^(-1)*s + w ),


% one less active constraint
Ma = [1 -3 -2; 2 1 2];
wa = [1;1];
lbd = -( Ma*2*Q^(-1)*Ma' )^(-1)*(4*Ma*Q^(-1)*s + wa ),

% final solution (one less active constraint)
Ma = [1 1 1];
wa = 1;
lbd_opt = -(2* Ma*Q^(-1)*Ma' )^(-1)*(4*Ma*Q^(-1)*s + wa ),
x_opt = -Q^-1*(Ma'*lbd_opt + s), 

% get numerial solution for comparison:
x_opt2 = quadprog(Q,s,M,w),

%% Exercise 2.6: Double integrator constrained MPC regulator

clear X U Xopt Uopt TX TU TXopt TUopt X2 U2;

A = [1 1;0 1];
B = [0;1];
N = 3;
P = eye(2);
Q = eye(2);
R = 4;
x0 = [-4.5;2];
nx = size(B,1);
nu = size(B,2);

% compute the batch matrices
[F,G,Qb,Rb,H] = GetBatchXMatrices(A,B,[],N,P,Q,R);
Qt = F'*Qb*F,
Rt = G'*Qb*G + Rb,
St = G'*Qb*F,
K = Rt^-1*St; % unconstrained gain (for comparison)

% compute constraints matrices:
u_max = 0.3;
x_max = [5;5];
U_max = kron(u_max,ones(N,1));
X_max = kron(x_max,ones(N+1,1));
Mu = [-eye(N*nu);eye(N*nu)];
wu = [U_max;U_max];
Mx = [-G;G];
M = [Mu;Mx];

%simulate controlled system:
X(:,1) = x0;
X2(:,1) = x0;
nk = 30;
TU = 1:nk;
TX = 1:nk+1;
for k = 1:nk
    % initializations:
    xk = X(:,k);
    wx = [X_max + F*xk;X_max-F*xk];
    w = [wu;wx];
    
    % compute constrained optimal sequence and MPC policy
    [Uo,Jo,exitflag,output,lambda] = quadprog(Rt,2*St*xk,M,w);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end    
    Uopt(:,:,k) = reshape( Uo ,nu,N);    
    U(:,k) = Uopt(:,1,k);
    
    % simulate system
    X(:,k+1) = A*X(:,k) + B*U(:,k);
    
    % compute unconstrained optimal sequence and MPC policy (for comparison):
    Uopt2(:,:,k) = reshape( -K*X2(:,k) ,nu,N);
    U2(:,k) = Uopt2(:,1,k);
    
    % simulate system 2 (for comparison):
    X2(:,k+1) = A*X2(:,k) + B*U2(:,k);
    
    % Compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( F*x0 + G*Uo ,nx,N+1);
    Xopt2(:,:,k) = reshape( (F-G*K)*X2(:,k) ,nx,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
end

figure(4301);
plot(X(1,:),X(2,:),'s-','Color',sstblue);
grid on;
hold on;
plot(X2(1,:),X2(2,:),'o-','Color',sstgray);
for k = 1:nk
    plot(Xopt(1,:,k),Xopt(2,:,k),'.--','Color',sstlightblue);
    plot(Xopt2(1,:,k),Xopt2(2,:,k),'.-.','Color',sstlightgray);
end
plot(X(1,:),X(2,:),'s-','Color',sstblue);
plot(X2(1,:),X2(2,:),'o-','Color',sstgray);
hold off;
xlabel('x_1');
ylabel('x_2');
legend('const.','unconst.','const. pred.','unconst. pred.');
title('Phase plot');

figure(4302);
plot(TX,X(1,:),'s-','Color',sstblue);
grid on;
hold on;
plot(TX,X(2,:),'d-','Color',sstgreen);
plot(TX,X2(1,:),'o-','Color',sstgray);
plot(TX,X2(2,:),'+-','Color',sstlightgray);
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-.','Color',sstlightgreen);
end
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
legend('x_1 const.','x_2 const.','x_1 unc.','x_2 unc.','x_1 pred.','x_2 pred.','Location','SouthEast');
title('State evolution');

figure(4303);
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
plot(TU,U2,'o-','Color',sstgray);
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.--','Color',sstlightblue);
    plot(TUopt(k,:),Uopt2(:,:,k),'.-.','Color',sstlightgray);
end
plot(TU,U,'s-','Color',sstblue);
plot(TU,U2,'o-','Color',sstgray);
hold off;
xlabel('k');
ylabel('u(k)');
legend('const.','unc.','const. pred.','unconst. pred.');
title('Input');

%% Exercise 2.7: Double integrator constrained tracking MPC
% Exercise 2.7

clear X U Xopt Uopt TX TU X2 U2 Xopt2 Uopt2;
example_name = 'doubleint_const_track_MPC_int_simul';

Ad = [1 1;0 1];
Bd = [0;1];
Cd = [1,0];
P = 1; %notice P has dimensions compatible with the output
Q = 1; %notice Q has dimensions compatible with the output
R = 4;
N = 3;
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

% compute constraints matrices:
u_max = 0.3;
y_max = 5;
U_max = kron(u_max,ones(N,1));
Y_max = kron(y_max,ones(N+1,1));
M3 = tril(ones(N*nu));
M4 = ones(N*nu,nu);
Mu = [-M3;M3];
My = [-Gb;Gb];
M = Mu;

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
    wu = [U_max + M4*u_1;U_max - M4*u_1];
    w = wu;
    Dxdk2 = Xd2(:,k)-Xd2(:,k-1);
    xk2 = [ Dxdk2; Cd*Xd2(:,k)];
    
    % compute constrained optimal incremental control sequence and MPC policy
    [dUo,Jo,exitflag,output,lambda] = quadprog(Rt,2*St*(Fb*xk-Yb),M,w);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end    
    dUopt(:,:,k) = reshape( dUo ,nu,N);    
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    U(:,k) = Uopt(:,1,k);
    
    % simulate original system:
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
    
figure(4401);
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

figure(4402);
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

figure(4403);
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

%% Exercise 2.8: no solution provided.
clear all;
sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;
clear X U Xopt Uopt TX TU X2 U2 Xopt2 Uopt2;
A = [1.3 1;0 1];
B = [0;1];
N = 3;
P = eye(2);
Q = eye(2);
R = 15;
x0 = [-2;2];
nx = size(B,1);
nu = size(B,2);

% compute the batch matrices
[F,G,Qb,Rb,H] = GetBatchXMatrices(A,B,[],N,P,Q,R);
Qt = F'*Qb*F,
Rt = G'*Qb*G + Rb,
St = G'*Qb*F,
K = Rt^-1*St; % unconstrained gain (for comparison)

% compute constraints matrices:
u_max = 0.8;
x_max = [6;6];
U_max = kron(u_max,ones(N,1));
X_max = kron(x_max,ones(N+1,1));
Mu = [-eye(N*nu);eye(N*nu)];
wu = [U_max;U_max];
Mx = [-G;G];
M = [Mu;Mx];

%simulate controlled system:
X(:,1) = x0;
X2(:,1) = x0;
nk = 30;
TU = 1:nk;
TX = 1:nk+1;
for k = 1:nk
    % initializations:
    xk = X(:,k);
    wx = [X_max + F*xk;X_max-F*xk];
    w = [wu;wx];
    
    % compute constrained optimal sequence and MPC policy
    [Uo,Jo,exitflag,output,lambda] = quadprog(2*Rt,2*St*xk,M,w);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end    
    Uopt(:,:,k) = reshape( Uo ,nu,N);    
    U(:,k) = Uopt(:,1,k);
    
    % simulate system
    X(:,k+1) = A*X(:,k) + B*U(:,k);
    
%     % compute unconstrained optimal sequence and MPC policy (for comparison):
%     Uopt2(:,:,k) = reshape( -K*X2(:,k) ,nu,N);
%     U2(:,k) = Uopt2(:,1,k);
%     
%     % simulate system 2 (for comparison):
%     X2(:,k+1) = A*X2(:,k) + B*U2(:,k);
    
    % Compute auxiliary variables for visualization:
    Xopt(:,:,k) = reshape( F*x0 + G*Uo ,nx,N+1);
    %Xopt2(:,:,k) = reshape( (F-G*K)*X2(:,k) ,nx,N+1);
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
end

figure(4301);
plot(X(1,:),X(2,:),'s-','Color',sstblue);
grid on;
hold on;
%plot(X2(1,:),X2(2,:),'o-','Color',sstgray);
for k = 1:nk
    plot(Xopt(1,:,k),Xopt(2,:,k),'.--','Color',sstlightblue);
    %plot(Xopt2(1,:,k),Xopt2(2,:,k),'.-.','Color',sstlightgray);
end
plot(X(1,:),X(2,:),'s-','Color',sstblue);
%plot(X2(1,:),X2(2,:),'o-','Color',sstgray);
hold off;
xlabel('x_1');
ylabel('x_2');
%legend('const.','unconst.','const. pred.','unconst. pred.');
legend('const.','const. pred.');
title('Phase plot');

figure(4302);
plot(TX,X(1,:),'s-','Color',sstblue);
grid on;
hold on;
plot(TX,X(2,:),'d-','Color',sstgreen);
%plot(TX,X2(1,:),'o-','Color',sstgray);
%plot(TX,X2(2,:),'+-','Color',sstlightgray);
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-.','Color',sstlightgreen);
end
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
hold off;
xlabel('k');
ylabel('x(k)');
%legend('x_1 const.','x_2 const.','x_1 unc.','x_2 unc.','x_1 pred.','x_2 pred.','Location','SouthEast');
legend('x_1 const.','x_2 const.','x_1 pred.','x_2 pred.','Location','SouthEast');
title('State evolution');

figure(4303);
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
%plot(TU,U2,'o-','Color',sstgray);
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.--','Color',sstlightblue);
    %plot(TUopt(k,:),Uopt2(:,:,k),'.-.','Color',sstlightgray);
end
plot(TU,U,'s-','Color',sstblue);
%plot(TU,U2,'o-','Color',sstgray);
hold off;
xlabel('k');
ylabel('u(k)');
%legend('const.','unc.','const. pred.','unconst. pred.');
legend('const.','const. pred.');
title('Input');
%% Exercise 2.9: no solution provided.


clear X U Xopt Uopt TX TU X2 U2 Xopt2 Uopt2;
example_name = 'doubleint_const_track_MPC_int_simul';

Ad = [0.3 0.6;0 1.2];
Bd = [0;1];
Cd = [1,0];
P = 1; %notice P has dimensions compatible with the output
Q = 1; %notice Q has dimensions compatible with the output
R = 4;
N = 3;
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

% compute constraints matrices:
u_max = 0.8;
y_max = 5;
U_max = kron(u_max,ones(N,1));
Y_max = kron(y_max,ones(N+1,1));
M3 = tril(ones(N*nu));
M4 = ones(N*nu,nu);
Mu = [-M3;M3];
My = [-Gb;Gb];
M = Mu;

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
    wu = [U_max + M4*u_1;U_max - M4*u_1];
    w = wu;
    %Dxdk2 = Xd2(:,k)-Xd2(:,k-1);
    %xk2 = [ Dxdk2; Cd*Xd2(:,k)];
    
    % compute constrained optimal incremental control sequence and MPC policy
    [dUo,Jo,exitflag,output,lambda] = quadprog(2*Rt,2*St*(Fb*xk-Yb),M,w);
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
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; Cd*Xd(:,k+1)];
    
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



%%

