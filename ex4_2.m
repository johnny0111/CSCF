clear all
close all
clc



Ad=[0.6 0;0 0.9];
Bd=[1.2,0.8;-0.8,1.2];

Ad_d=[0.6 0;0 0.9];
Bd_d=[1.2 0;0 1.2]

Cd=[1 0;0 1];
N=3;
x_inicial=[1;1];
Ts=0.1;

y_ref_1=[repmat(0,1,30) repmat(1,1,30)];
y_ref_2=[repmat(0,1,30) repmat(1.5,1,30)];

R=eye(2,2)*1;
P=eye(2,2)*5;
Q=eye(2,2)*5;


[F,G,Qb,Rb,H,Fd,Gd,Hd,A,B,C] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);


F_l=H*F;
G_l=H*G;

R_solu=G_l'*Qb*G_l+Rb;
S_solu=G_l'*Qb;
K_y=inv(R_solu)*S_solu;
K_solu=K_y*F_l;

[F_d,G_d,Qb_d,Rb_d,H_d,Fd_d,Gd_d,Hd_d,A_d,B_d,C_d] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
F_l_d=H_d*F;
G_l_d=H_d*G;

R_solu_d=G_l'*Qb*G_l+Rb;
S_solu_d=G_l'*Qb;
K_y_d=inv(R_solu_d)*S_solu_d;
K_solu_d=K_y_d*F_l_d;

x0=[x_inicial*0 ; Cd*x_inicial];
nu=2;
nk=60;

U = zeros(nu,nk);
dU = [];

U_d = zeros(nu,nk);
dU_d = [];

Xd(:,1) = x_inicial;
X(:,1) = x0;
Y(:,1) = Cd*x_inicial;
Xd(:,2) = x_inicial;
X(:,2) = x0;
Y(:,2) = Cd*x_inicial;


Xd_d(:,1) = x_inicial;
Xd_d(:,2) = x_inicial;
X_d(:,1)=x0;
X_d(:,2)=x0;
Y_d(:,1) = Cd*x_inicial;
Y_d(:,2) = Cd*x_inicial;

for k=2:57
    
    %Define first conditions
    %Yb=[y_ref_1(:,k:k+N) y_ref_2(:,k:k+N)]';
    Yb=[y_ref_1(:,k+0) y_ref_2(:,k+0) y_ref_1(:,k+1) y_ref_2(:,k+1) y_ref_1(:,k+2) y_ref_2(:,k+2) y_ref_1(:,k+3) y_ref_2(:,k+3)]'
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; Cd*Xd(:,k)];
    
    Yb_d=[y_ref_1(:,k+0) y_ref_2(:,k+0) y_ref_1(:,k+1) y_ref_2(:,k+1) y_ref_1(:,k+2) y_ref_2(:,k+2) y_ref_1(:,k+3) y_ref_2(:,k+3)]'
    Dxdk_d = Xd_d(:,k)-Xd_d(:,k-1);
    X_d(:,k) = [ Dxdk_d; Cd*Xd_d(:,k)];
    
    %Controlador
    dU_aux=-K_solu*X(:,k)+K_y*Yb;
    dU=[dU dU_aux]
    
    dU_aux_d=-K_solu_d*X_d(:,k)+K_y_d*Yb;
    dU_d=[dU_d dU_aux_d]
    
    %Acao de controlo
    dU0=[dU_aux(1,1);dU_aux(2,1)];
    U(:,k)=dU0+U(:,k-1);
    
    dU0_d=[dU_aux_d(1,1);dU_aux_d(2,1)];
    U_d(:,k)=dU0_d+U_d(:,k-1);

    %Instalacao
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k) ;
    Y(:,k+1) = Cd*Xd(:,k+1);
    
    Xd_d(:,k+1) = Ad_d*Xd_d(:,k) + Bd_d*U_d(:,k) ;
    Y_d(:,k+1) = Cd*Xd_d(:,k+1);
   
end    

figure;
title("Centralized Regulation")
t=Ts:0.1:6
plot(t(3:length(t)),Y)
hold on
plot(t,[y_ref_1' y_ref_2'])
legend('Player 1 Output','Player 2 Output','Player 1 reference','Player 2 reference')


figure;
title("Decentralized Regulation")
t=Ts:0.1:6
plot(t(3:length(t)),Y_d)
hold on
plot(t,[y_ref_1' y_ref_2'])
legend('Player 1 Output','Player 2 Output','Player 1 reference','Player 2 reference')