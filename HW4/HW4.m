%%
clc;
clear all;
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

%% Ex 4.3 com efeito integral e sem restrições

np = 2;
% example_name = 'exercise43_distributed_doubleint_p10';
% np = 10;

Ai = [  1 , 0.1 , -0.1 
        0 ,   1 ,  0
		0 ,   0 ,  1 ] ;
A1 = Ai;
A2 = Ai;
A = Ai;
B11 = [ 0 ; 1 ; 0];
B22 = [ 0 ; 0 ; 1];
B12 = B22;
B21 = B11;
B = [B11,B22];
C1 = [ 1 0 0];
C2 = [-1 0 0];
C = [C1;C2];

x10 = [1;1;2];
x20 = x10;
x0 = x10;

nu = size(B,2);
ny = size(C,1);
nx = size(B,1);
nx1 = size(B11,1);
nx2 = size(B22,1);
nu1 = size(B11,2);
nu2 = size(B22,2);
ny1 = size(C1,1);
ny2 = size(C2,1);

A1_ext = [A1 zeros(nx1,ny1);C1*A1 eye(ny1)];
A2_ext = [A2 zeros(nx1,ny1);C2*A2 eye(ny2)];
A_ext = [A zeros(nx,ny);C*A eye(ny)];
B11_ext = [ B11;C1*B11];
B22_ext = [ B22;C2*B22];
B12_ext = [ B12;C1*B12];
B21_ext = [ B21;C2*B21];
B_ext = [B;C*B];
C1_ext = [zeros(ny1,nx1),eye(ny1)];
C2_ext = [zeros(ny2,nx2),eye(ny2)];
C_ext = [zeros(ny,nx),eye(ny)];




% cost parameters
N = 3;
Pi = 2;
Qi = 2;
Ri = 1;
alpha1 = 1;
alpha2 = 1;

% distributed steps parameters
w1 = 0.4;
w2 = 1-w1;

% compute centralized tracking controller
P = blkdiag(Pi,Pi);
Q = blkdiag(Qi,Qi);
R = blkdiag(Ri,Ri);
[Fb,Gb,Qb,Rb,F,G,H] = GetBatchYMatrices(A_ext,B_ext,C_ext,N,P,Q,R);
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Kcy = Rt^-1*St;
Kc = -Kcy*Fb;


% compute decentralized tracking controllers
[Fb1,Gb1,Qb1,Rb1,F1,G1,H1] = GetBatchYMatrices(A1_ext,B11_ext,C1_ext,N,Pi,Qi,Ri);
Rt1 = Gb1'*Qb1*Gb1 + Rb1;
St1 = Gb1'*Qb1;
Kdy1 = Rt1^(-1)*St1;
Kd1 = -Kdy1*Fb1;

[Fb2,Gb2,Qb2,Rb2,F2,G2,H2] = GetBatchYMatrices(A2_ext,B22_ext,C2_ext,N,Pi,Qi,Ri);
Rt2 = Gb2'*Qb2*Gb2 + Rb2;
St2 = Gb2'*Qb2;
Kdy2 = Rt2^(-1)*St2;
Kd2 = -Kdy2*Fb2;
% compute decentralized tracking controllers
% [Fb1,Gb1,Qb1,Rb1,F1,G1,H1] = GetBatchYMatrices(A1,B11,C1,N,Pi,Qi,Ri);
% Rt1 = Gb1'*Qb1*Gb1 + Rb1;
% St1 = Gb1'*Qb1;
% Kdy1 = Rt1^(-1)*St1;
% Kd1 = -Kdy1*Fb1;
% 
% [Fb2,Gb2,Qb2,Rb2,F2,G2,H2] = GetBatchYMatrices(A2,B22,C2,N,Pi,Qi,Ri);
% Rt2 = Gb2'*Qb2*Gb2 + Rb2;
% St2 = Gb2'*Qb2;
% Kdy2 = Rt2^(-1)*St2;
%  Kd2 = -Kdy2*Fb2;

[Fbd,Gbd,Qbd,Rbd,Fd,Gd,Hd] = GetYMats(A_ext,B_ext,C_ext,N,P,Q,R,alpha1);
[Rt,Sx,Sy,Su,K,Ky,L] = GetYMPC(Fbd,Gbd,Qbd,Rbd);

nk = 100;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = [2*(Tref>=50);-2*(Tref>=50)]; % step reference of amplitude 1, starting at k=30
uc = zeros(nu,nk);
ud = zeros(nu,nk);
u = zeros(nu,nk);
U1 = zeros(nu1,N,nk);
U2 = zeros(nu2,N,nk);

% x_ = [x0*0 ; C*x0];
xc(:,1) = x0;
% X(:,1) = x_;
xc(:,2) = x0;


xd(:,1) = x0;
xd(:,2) = x0;
% xd2(:,1) = x20;
% xd2(:,2) = x20;
for k = 2:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    
    % centralized MPC
    Dxdk = xc(:,k)-xc(:,k-1);
    X(:,k) = [ Dxdk; C*xc(:,k)];
    xk = X(:,k);
    Ucopt(:,:,k) = reshape( Kc*xk+Kcy*Yb ,nu,N);
    u_c= Ucopt(:,1,k);
    uc(:,k) =u_c + Ucopt(:,1,k-1);

    xc(:,k+1) = A*xc(:,k) + B*uc(:,k); %simulate joint system
    
    
   % decentralized MPC
    Dxdkd1 = xd(:,k)-xd(:,k-1);
    Xd(:,k) = [Dxdkd1; C*xd(:,k)];
    

    
    Ud1opt(:,:,k) = reshape( Kd1*[Xd(1:3,k);Xd(4,k)]+Kdy1*Yb1 ,nu1,N);
    Ud2opt(:,:,k) = reshape( Kd2*[Xd(1:3,k);Xd(5,k)]+Kdy2*Yb2 ,nu2,N);
    ud1 = Ud1opt(:,1,k);    
    ud2 = Ud2opt(:,1,k);    
    Dud(:,k) = [ ud1 ; ud2 ];
    ud(:,k) = Dud(:,k) + ud(:,k-1);
    xd(:,k+1) = A*xd(:,k) + B*ud(:,k); %simulate joint system
%        xd1 = xd(:,k);
%     xd2 = xd(:,k);
%     Ud1opt(:,:,k) = reshape( Kd1*xd1+Kdy1*Yb1 ,nu1,N);
%     Ud2opt(:,:,k) = reshape( Kd2*xd2+Kdy2*Yb2 ,nu2,N);
%     ud1 = Ud1opt(:,1,k);    
%     ud2 = Ud2opt(:,1,k);    
%     ud(:,k) = [ ud1 ; ud2 ];
%     
%     % simulate system for decentralized MPC
%     xd(:,k+1) = A*xd(:,k) + B*ud(:,k); %simulate joint system
    
    
end



figure(7201);
plot(ref(1,:),ref(2,:),'k+-');
grid on;
hold on;
plot(xc(1,:),xc(2,:),'s-','Color',sstblue);
plot(xc(1,:),xc(3,:),'.-','Color',sstdarkblue);
plot(xd(1,:),xd(2,:),'o-','Color',sstgray);
plot(xd(1,:),xd(3,:),'.-','Color',sstlightgray);
% plot(x(1,:),x(2,:),'x-','Color',sstgreen);
% plot(x(1,:),x(3,:),'.--','Color',sstdarkgreen);
hold off;
xlabel('$$\theta_{12}$$');
ylabel('$$\omega_i$$');
legend('ref.','cent. 1','cent. 2','decent. 1','decent. 2','dist. 1','dist. 2');
title('Phase plot');

figure(7202);
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
% plot(TX,x(1,:),'o-','Color',sstdarkgreen);
% plot(TX,x(2,:),'.--','Color',sstgreen);
% plot(TX,x(3,:),'.-.','Color',sstlightgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r_1$$','$$r_2$$','$$x_1$$ cen.','$$x_2$$ cen.','$$x_3$$ cen.','$$x_1$$ decen.','$$x_2$$ decen.','$$x_3$$ decen.','$$x_1$$ dist.','$$x_2$$ dist.','$$x_3$$ dist.');
title('State evolution');

figure(7203);
plot(TU,uc(1,:),'s-','Color',sstdarkblue);
grid on;
hold on;
plot(TU,uc(2,:),'.--','Color',sstblue);
% plot(TU,ud(1,:),'d-','Color',sstdarkgray);
% plot(TU,ud(2,:),'.-.','Color',sstgray);
% plot(TU,u(1,:),'o-','Color',sstdarkgreen);
% plot(TU,u(2,:),'.--','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$ cent.','$$u_2$$ cent.','$$u_1$$ decent.','$$u_2$$ decent.','$$u_1$$ dist.','$$u_2$$ dist.');
title('Input');















