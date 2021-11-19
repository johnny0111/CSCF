clear all;
sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;
clear X U Xopt Uopt TX TU X2 U2 Xopt2 Uopt2;
A = [1.5 1;0 0.4];
B = [0;1];
N = 3;
P = eye(2);
Q = eye(2);
R = 8;
x0 = [-3;1];
nx = size(B,1);
nu = size(B,2);

% compute the batch matrices
[F,G,Qb,Rb,H] = GetBatchXMatrices(A,B,[],N,P,Q,R);
Qt = F'*Qb*F,
Rt = G'*Qb*G + Rb,
St = G'*Qb*F,
K = Rt^-1*St; % unconstrained gain (for comparison)

% compute constraints matrices:
u_min = -0.5;
u_max = 1.5;
x_max = [4;4];
U_max = kron(u_max,ones(N,1));
U_min = kron(u_min,ones(N,1));
X_max = kron(x_max,ones(N+1,1));
Mu = [-eye(N*nu);eye(N*nu)];
wu = [-U_min;U_max];
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

%%
Q = eye(3);
s = [-1;-2;1];
M = [1 -3 -2; 2 1 2; 1 -2 2];
w = [-1;1;2];

% check for inactive constraints:
det(M),
lbd = -( M*2*Q^(-1)*M' )^(-1)*(M*4*Q^(-1)*s + w ),

% one less active constraint
Ma = [1 -3 -2; 2 1 2];
wa = [-1;1];
lbd = -( M*2*Q^(-1)*M' )^(-1)*(M*4*Q^(-1)*s + w ),
return
% final solution (one less active constraint)
Ma = [2 1 2];
wa = 1;
lbd_opt = -( M*2*Q^(-1)*M' )^(-1)*(M*4*Q^(-1)*s + w),
x_opt = -Q^-1*(Ma'*lbd_opt + s), 

% get numerial solution for comparison:
x_opt2 = quadprog(Q,s,M,w),