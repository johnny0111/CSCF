%%
clearvars ;
% close all;
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
%% setup & model
p=5;
m = 1500;
Cd = 0.3;
area = 1.8;
rho = 1.225;
g = 9.8065;
theta = 0.005;
fmax = 2800;
fmin = -2200;
deltaFmax = 200;%200;%300;
maxTheta = 0.1;
dDist = 10;
distMin = 3;
vRe = 25;
ve = 25;
ws = 15;
Ts = 0.1;
d = (rho*area*Cd*ve)/m;
vl = 0;
lambda = 0.95;

% Ac =[0  -1  0   0   0   0   0   0   0   0   0   0   0   0   
%      0  -d  0   0   0   0   0   0   0   0   0   0   0   0   
%      0  1   0   -1  0   0   0   0   0   0   0   0   0   0   
%      0  0   0   -d  0   0   0   0   0   0   0   0   0   0   
%      0  0   0   1   0   -1  0   0   0   0   0   0   0   0   
%      0  0   0   0   0   -d  0   0   0   0   0   0   0   0   
%      0  0   0   0   0   1   0   -1  0   0   0   0   0   0   
%      0  0   0   0   0   0   0   -d  0   0   0   0   0   0   
%      0  0   0   0   0   0   0   1   0   -1  0   0   0   0   
%      0  0   0   0   0   0   0   0   0   -d  0   0   0   0   
%      0  0   0   0   0   0   0   0   0   1   0   -1  0   0  
%      0  0   0   0   0   0   0   0   0   0   0   -d  0   0   
%      0  0   0   0   0   0   0   0   0   0   0   1   0   -1 
%      0  0   0   0   0   0   0   0   0   0   0   0   0   -d  ];
%  Ac = [0  -1  0   0       
%        0  -d  0   0       
%        0  1   0   -1     
%        0  0   0   -d    ];  
% 
%  Ac = [0  -1  0   0     0   0
%        0  -d  0   0     0   0   
%        0  1   0   -1    0   0  
%        0  0   0   -d    0   0      
%        0  0   0    1    0   -1
%        0  0   0    0    0   -d];    
Ac =[0  -1  0   0   0   0   0   0   0   0   
     0  -d  0   0   0   0   0   0   0   0   
     0  1   0   -1  0   0   0   0   0   0   
     0  0   0   -d  0   0   0   0   0   0  
     0  0   0   1   0   -1  0   0   0   0   
     0  0   0   0   0   -d  0   0   0   0   
     0  0   0   0   0   1   0   -1  0   0  
     0  0   0   0   0   0   0   -d  0   0   
     0  0   0   0   0   0   0   1   0   -1    
     0  0   0   0   0   0   0   0   0   -d ];
 
 Bd = [1  0  0
       0  -g  -d/2
       0   0  0
       0  -g  -d/2
       0   0  0
       0  -g  -d/2
       0   0  0
       0  -g  -d/2
       0   0  0
       0  -g  -d/2];
   

% Bc= [  0   0    0   0   0   0   0   
%        1/m 0    0   0   0   0   0   
%         0   0   0   0   0   0   0   
%         0  1/m  0   0   0   0   0   
%         0   0   0   0   0   0   0  
%         0   0   1/m 0   0   0   0   
%         0   0   0   0   0   0   0   
%         0   0   0   1/m 0   0   0   
%         0   0   0   0   0   0   0   
%         0   0   0   0   1/m 0   0   
%         0   0   0   0   0   0   0   
%         0   0   0   0   0   1/m 0   
%         0   0   0   0   0   0   0   
%         0   0   0   0   0   0   1/m ];
%     
%     Bc= [  0   0       
%            1/m 0        
%             0   0       
%             0  1/m   ];
% Bc= [  0   0   0     
%        1/m 0   0       
%         0   0  0      
%         0  1/m 0      
%         0   0  0      
%         0   0  1/m 
%  ];
Bc= [  0   0    0   0   0   
       1/m 0    0   0   0   
        0   0   0   0   0   
        0  1/m  0   0   0   
        0   0   0   0   0   
        0   0   1/m 0   0   
        0   0   0   0   0   
        0   0   0   1/m 0    
        0   0   0   0   0    
        0   0   0   0   1/m  ];
% Bijc = {Bc(1:2,1),Bc(1:2,2),Bc(1:2,3),Bc(1:2,4),Bc(1:2,5),Bc(1:2,6),Bc(1:2,7)
%         Bc(2:4,1),Bc(2:4,2),Bc(2:4,3),Bc(2:4,4),Bc(2:4,5),Bc(2:4,6),Bc(2:4,7)
%         Bc(4:6,1),Bc(4:6,2),Bc(4:6,3),Bc(4:6,4),Bc(4:6,5),Bc(4:6,6),Bc(4:6,7)
%         Bc(6:8,1),Bc(6:8,2),Bc(6:8,3),Bc(6:8,4),Bc(6:8,5),Bc(6:8,6),Bc(6:8,7)
%         Bc(8:10,1),Bc(8:10,2),Bc(8:10,3),Bc(8:10,4),Bc(8:10,5),Bc(8:10,6),Bc(8:10,7)
%         Bc(10:12,1),Bc(10:12,2),Bc(10:12,3),Bc(10:12,4),Bc(10:12,5),Bc(10:12,6),Bc(10:12,7)
%         Bc(12:14,1),Bc(12:14,2),Bc(12:14,3),Bc(12:14,4),Bc(12:14,5),Bc(12:14,6),Bc(12:14,7)};
%     Bijc = {Bc(1:2,1),Bc(1:2,2)
%         Bc(2:4,1),Bc(2:4,2)};

% Bijc = {Bc(1:2,1),Bc(1:2,2),Bc(1:2,3)
%         Bc(2:4,1),Bc(2:4,2),Bc(2:4,3)
%         Bc(4:6,1),Bc(4:6,2),Bc(4:6,3)
%        };

Bijc = {Bc(1:2,1),Bc(1:2,2),Bc(1:2,3),Bc(1:2,4),Bc(1:2,5)
        Bc(2:4,1),Bc(2:4,2),Bc(2:4,3),Bc(2:4,4),Bc(2:4,5)
        Bc(4:6,1),Bc(4:6,2),Bc(4:6,3),Bc(4:6,4),Bc(4:6,5)
        Bc(6:8,1),Bc(6:8,2),Bc(6:8,3),Bc(6:8,4),Bc(6:8,5)
        Bc(8:10,1),Bc(8:10,2),Bc(8:10,3),Bc(8:10,4),Bc(8:10,5)};
% C=[ 1  0   0   0   0   0   0   0   0   0   0   0   0   0   
%     0  0   1   0   0   0   0   0   0   0   0   0   0   0   
%     0  0   0   0   1   0   0   0   0   0   0   0   0   0   
%     0  0   0   0   0   0   1   0   0   0   0   0   0   0   
%     0  0   0   0   0   0   0   0   1   0   0   0   0   0  
%     0  0   0   0   0   0   0   0   0   0   1   0   0   0   
%     0  0   0   0   0   0   0   0   0   0   0   0   1   0   ];
% C=[ 1  0   0   0      
%     0  0   1   0  ];
% C=[ 1  0   0   0   0   0   
%     0  0   1   0   0   0   
%     0  0   0   0   1   0];

C=[ 1  0   0   0   0   0   0   0   0   0  
    0  0   1   0   0   0   0   0   0   0  
    0  0   0   0   1   0   0   0   0   0   
    0  0   0   0   0   0   1   0   0   0   
    0  0   0   0   0   0   0   0   1   0  ];

% Ci={[1 0], [0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0]};
% Ci={[1 0], [0 1 0]};
%  Ci={[1 0], [0 1 0],[0 1 0]};
Ci={[1 0], [0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0]};

 A1c = [0 -1;0 -d];
 A2c = [-d 0 0;1 0 -1;0 0 -d];
 
% A = eye(6) + Ac*Ts;
% Ai={eye(2) + A1c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts,...
%      eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts};
%  Ai={eye(2) + A1c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts};
Ai={eye(2) + A1c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts};
A=eye(2*p) + Ac*Ts;
% Ai={eye(2) + A1c*Ts, eye(3) + A2c*Ts};
B = Bc*Ts;
for i=1:p
    for j=1:p
       Bij{i,j} = Bijc{i,j}*Ts;
    end
    
end


nu = size(B,2);
ny = size(C,1);
nx = size(B,1);
nx1 = size(Bij{1,1},1);
nu1 = size(Bij{1,1},2);
nx2 = size(Bij{2,1},1);
nu2 = size(Bij{2,1},2);
ny1 = size(Ci{1},1);
ny2 = size(Ci{2},1);
% x10 = [1;1];
% x20 = [1;1;1];
% x30 = [1;1;1];
% x40 = [1;1;1];
% x50 = [1;1;1];
% x60 = [1;1;1];
% x70 = [1;1;1];
% x80 = [1;1;1];
% x90 = [1;1;1];
% x100 = [1;1;1];

% xd0=[x10;x20;x30;x40;x50;x60;x70;x80;x90;x100];
xd0 = ones(2*p,1);


%A_ext = [A zeros(nx,ny);C*A eye(ny)];
Ai_ext{1} = [Ai{1} zeros(nx1,ny1);Ci{1}*Ai{1} eye(ny1)];
for i=2:p
   Ai_ext{i} = [Ai{i} zeros(nx2,ny2);Ci{i}*Ai{i} eye(ny2)];
end

for i=1:p
   Bij_ext{1,i} = [Bij{1,i};Ci{1}*Bij{1,i}]; 
end
 
for i=2:p
   for j = 1:p
        Bij_ext{i,j} = [Bij{i,j};Ci{i}*Bij{i,1}]; 
   end
end
 
Ci_ext{1} = [zeros(ny1,nx1),eye(ny1)];

for i=2:p
   Ci_ext{i} = [zeros(ny2,nx2),eye(ny2)]; 
end
 
 %% Cost
 N = 60;%60 p =3 ;
Pi = 9999999999;
Pi = 999;
% qi=50000; %p=2
% ri=10; %p=2
% 
% Qi = [qi; 0.8*qi]; 
% Ri = [ri; 0.8*ri];

% qi=0.1;%p=3
% ri = 0.00001;%p=3

% Qi = [qi; 0.5*qi; 0.005*qi];
% Ri = [ri; 0.1*ri; 0.03*ri];


% qi=0.1;%p=4
% ri = 0.00001;%p=4
% 
% Qi = [qi; 0.5*qi; 0.005*qi; 0.0005*qi];
% Ri = [ri; 0.1*ri; 0.06*ri; 0.003*ri];

qi=0.1;%p=5
ri = 0.00001;%p=5

Qi = [qi; 0.5*qi; 0.005*qi; 0.0006*qi; 0.00005*qi];
Ri = [ri; 0.1*ri; 0.06*ri; 0.005*ri; 0.0000005*ri];

% qi=0.1;%p=6
% ri = 0.00001;%p=6
% 
% Qi = [qi; 0.5*qi; 0.005*qi; 0.0006*qi; 0.00005*qi;  0.00001*qi];
% Ri = [ri; 0.1*ri; 0.06*ri; 0.005*ri; 0.0000005*ri; 0.0000001*ri];

alphai = [1 1 1 1 1 1 1 1 1 1 ];
% distributed steps parameters
w(1) = 1/p;
w(2) = 1/p;
w(3) = 1/p;
w(4) = 1/p;
w(5) = 1/p;
w(6) = 1/p;
w(7) = 1/p;
% w(8) = 1/10;
% w(9) = 1/10;
% w(10) = 1/10;
% for i=2:3
%    w(i) = (1-w(1))/9; 
% end
np=2;

%% Batch matrices
[Fb,Gb,Qb,Rb] = GetBatchYNdistMatrices(Ai_ext,Bij_ext,Ci_ext,N,Pi,Qi,Ri,alphai(1));
[Rt,Sx,Sy,Su,K,Ky,L] = GetYMPC(Fb,Gb,Qb,Rb);

% %% Constraints 
% 
%player 1

u_max = fmax -(0.5*rho*area*Cd*ve^2);
u_min = fmin - (0.5*rho*area*Cd*ve^2);
% U_max1 = kron(u_max1,ones(N*nu1,1));
% U_min1 = kron(u_min1,ones(N*nu1,1));
% M31 = kron(tril(ones(N)), eye(nu1));
% M41 = kron(ones(N,1), eye(nu1));
% Mu1 = [-M31;M31];
% 
du_max = deltaFmax;
du_min = -deltaFmax;
% DU_max1 = kron(du_max1,ones(N*nu1,1));
% DU_min1 = kron(du_min1,ones(N*nu1,1));
% Mdu1 = [-eye(N*nu1);eye(N*nu1)];
% wdu1 = [-DU_min1;DU_max1];

% pr_min = -7 ;
% pr_max = 93;
% Y_min1 = kron(pr_min,ones(nu1*(N+1),1));
% Y_max1 = kron(pr_max, ones(nu1*(N+1),1));
% My1 = [-Gb11; Gb11];
U_max = kron(u_max,ones(N,1));
U_min = kron(u_min,ones(N,1));
DU_max = kron(du_max,ones(N,1));
DU_min = kron(du_min,ones(N,1));
M3 = tril(ones(N*nu1));
M4 = ones(N*nu1,nu1);
Mu = [-M3;M3];
Mdu = [-eye(N);eye(N)];
wdu = [-DU_min;DU_max];
M1 = [Mdu; Mu];
% 
% 
% 
% 
% 
% %player i
% 
% u_maxi = fmax -(0.5*rho*area*Cd*ve^2);
% u_mini = fmin - (0.5*rho*area*Cd*ve^2);
% U_maxi = kron(u_maxi,ones(N*nu2,1));
% U_mini = kron(u_mini,ones(N*nu2,1));
% M3i = kron(tril(ones(N)), eye(nu2));
% M4i = kron(ones(N,1), eye(nu2));
% Mui = [-M3i;M3i];
% 
% du_maxi = deltaFmax;
% du_mini = -deltaFmax;
% DU_maxi = kron(du_maxi,ones(N*nu2,1));
% DU_mini = kron(du_mini,ones(N*nu2,1));
% Mdui = [-eye(N*nu2);eye(N*nu2)];
% wdui = [-DU_mini;DU_maxi];

% pr_min = -7 ;
% pr_max = 93;
% Y_min2 = kron(pr_min,ones(nu2*(N+1),1));
% Y_max2 = kron(pr_max, ones(nu2*(N+1),1));
% My2 = [-Gb22; Gb22];

U_maxi = kron(u_max,ones(N,1));
U_mini = kron(u_min,ones(N,1));
DU_maxi = kron(du_max,ones(N,1));
DU_mini = kron(du_min,ones(N,1));
M3i = tril(ones(N*nu1));
M4i = ones(N*nu2,nu2);
Mui = [-M3i;M3i];
Mdui = [-eye(N);eye(N)];
wdui = [-DU_mini;DU_maxi];

Mi = [Mdui; Mui];






%% Simulation

% close all;
nk = 800;%250;
TU = 1:nk;
TX = 1:nk+1;
%TX = 0:T:(nk+1)*T - T;
Tref = 1:nk+N;
ref = 10 * square(0.0002*Tref, 0.79);
%ref = [ref; ref];

dist = [zeros(1,nk+N);theta*(Tref>=500);ws*(Tref>=200)];

%player1
x01 = [xd0(1:2)*0 ; Ci{1}*xd0(1:2)];
U1 = zeros(nu1,N,nk);
X1(:,1) = x01;

%player2
x02 = [xd0(2:4)*0 ; Ci{2}*xd0(2:4)];
U2 = zeros(nu2,N,nk);
X2(:,1) = x02;

%player3
x03 = [xd0(4:6)*0 ; Ci{3}*xd0(4:6)];
U3 = zeros(nu2,N,nk);
X3(:,1) = x03;

%player4
x04 = [xd0(6:8)*0 ; Ci{4}*xd0(6:8)];
U4 = zeros(nu2,N,nk);
X4(:,1) = x04;

%player5
x05 = [xd0(8:10)*0 ; Ci{5}*xd0(8:10)];
U5 = zeros(nu2,N,nk);
X5(:,1) = x05;

% % player6
% x06 = [xd0(10:12)*0 ; Ci{6}*xd0(10:12)];
% U6 = zeros(nu2,N,nk);
% X6(:,1) = x06;

% %player7
% x07 = [xd0(12:14)*0 ; Ci{7}*xd0(12:14)];
% U7 = zeros(nu2,N,nk);
% X7(:,1) = x07;

% %player8
% x08 = [xd0(14:16)*0 ; Ci{8}*xd0(14:16)];
% U8 = zeros(nu2,N,nk);
% X8(:,1) = x08;
% 
% %player9
% x09 = [xd0(16:18)*0 ; Ci{9}*xd0(16:18)];
% U9 = zeros(nu2,N,nk);
% X9(:,1) = x09;
% 
% %player10
% x010 = [xd0(18:20)*0 ; Ci{10}*xd0(18:20)];
% U10 = zeros(nu2,N,nk);
% X10(:,1) = x010;

Xd(:,1) = xd0;
Xd(:,2) = xd0;
pk = eye(10);
u = zeros(p,nk);
Y = zeros(p,nk);
noise = zeros(p,nk);
c = 0.8;

for k=1:nk
    noise(:,k)=[c*rand;c*rand; c*rand; c*rand; c*rand; ]; 
end

for k=1:p
   V(k,1) = var(noise(k,:)); 
end
V=0.0000000001*V;

% sys = ss(A,B,C,0);
% sys = setDelayModel(sys,Ts*1.2);
% [A,B1,B2,C1,C2,D11,D12,D21,D22,E,tau] = getDelayModel(sys);
for k = 2:nk
    k
    Yb = ref(:,k:k+N)';

    %Player 1
    Dxdk1 = Xd(1:2,k)-Xd(1:2,k-1);
    X1(:,k) = [ Dxdk1; Ci{1}*Xd(1:2,k)];
    x1k = X1(:,k);
    
    %Player 2
    Dxdk2 = Xd(2:4,k)-Xd(2:4,k-1);
    X2(:,k) = [ Dxdk2; Ci{2}*Xd(2:4,k)];
    x2k = X2(:,k);
    
    %Player 3
    Dxdk3 = Xd(4:6,k)-Xd(4:6,k-1);
    X3(:,k) = [ Dxdk3; Ci{3}*Xd(4:6,k)];
    x3k = X3(:,k);
    
    %Player 4
    Dxdk4 = [Xd(6,k);Xd(7,k);Xd(8,k)]-[Xd(6,k-1);Xd(7,k-1);Xd(8,k-1)];
    X4(:,k) = [ Dxdk4; Ci{4}*[Xd(6,k);Xd(7,k);Xd(8,k)]];
    x4k = X4(:,k);
    
    %Player 5
    Dxdk5 = [Xd(8,k);Xd(9,k);Xd(10,k)]-[Xd(8,k-1);Xd(9,k-1);Xd(10,k-1)];
    X5(:,k) = [ Dxdk5; Ci{5}*[Xd(8,k);Xd(9,k);Xd(10,k)]];
    x5k = X5(:,k);    
    
%     %Player 6
%     Dxdk6 = [Xd(10,k);Xd(11,k);Xd(12,k)]-[Xd(10,k-1);Xd(11,k-1);Xd(12,k-1)];
%     X6(:,k) = [ Dxdk6; Ci{6}*[Xd(10,k);Xd(11,k);Xd(12,k)]];
%     x6k = X6(:,k);     
    
%     %Player 7
%     Dxdk7 = [Xd(12,k);Xd(13,k);Xd(14,k)]-[Xd(12,k-1);Xd(13,k-1);Xd(14,k-1)];
%     X7(:,k) = [ Dxdk7; Ci{7}*[Xd(12,k);Xd(13,k);Xd(14,k)]];
%     x7k = X7(:,k);    
    
%     %Player 8
%     Dxdk8 = [Xd(14,k);Xd(15,k);Xd(16,k)]-[Xd(14,k-1);Xd(15,k-1);Xd(16,k-1)];
%     X8(:,k) = [ Dxdk8; Ci{8}*[Xd(14,k);Xd(15,k);Xd(16,k)]];
%     x8k = X8(:,k); 
%     
%     %Player 9
%     Dxdk9 = [Xd(16,k);Xd(17,k);Xd(18,k)]-[Xd(16,k-1);Xd(17,k-1);Xd(18,k-1)];
%     X9(:,k) = [ Dxdk9; Ci{9}*[Xd(16,k);Xd(17,k);Xd(18,k)]];
%     x9k = X9(:,k);     
%     
%     %Player 10
%     Dxdk10 = [Xd(18,k);Xd(19,k);Xd(20,k)]-[Xd(18,k-1);Xd(19,k-1);Xd(20,k-1)];
%     X10(:,k) = [ Dxdk10; Ci{10}*[Xd(18,k);Xd(19,k);Xd(20,k)]];
%     x10k = X10(:,k);     
    
    wu1 = [-U_min + M4*u(1,k-1);U_max - M4*u(1,k-1)];
    wu2 = [-U_mini + M4i*u(2,k-1);U_maxi - M4i*u(2,k-1)];
    wu3 = [-U_mini + M4i*u(3,k-1);U_maxi - M4i*u(3,k-1)];
    wu4 = [-U_mini + M4i*u(4,k-1);U_maxi - M4i*u(4,k-1)];
    wu5 = [-U_mini + M4i*u(5,k-1);U_maxi - M4i*u(5,k-1)];
%     wu6 = [-U_mini + M4i*u(6,k-1);U_maxi - M4i*u(6,k-1)];
%     wu7 = [-U_mini + M4i*u(7,k-1);U_maxi - M4i*u(7,k-1)];
    %wu1 = [-U_min1 + M41*u(1,k-1);U_max1 - M41*u(1,k-1)];
    %wu2 = [-U_min2 + M42*u(2,k-1);U_max2 - M41*u(2,k-1)];
    %wy1 = [-Y_min1 + Fb1*x1k; Y_max1 - Fb1*x1k];
    %wy2 = [-Y_min2 + Fb2*x2k; Y_max2 - Fb2*x2k];
    wr1 = [wdu;wu1];
    wr2 = [wdui;wu2];
    wr3 = [wdui;wu3];
    wr4 = [wdui;wu4];
    wr5 = [wdui;wu5];
%     wr6 = [wdui;wu6];
%     wr7 = [wdui;wu7];
    
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    U3p = reshape( U3(:,:,k-1) ,[],1);
    U4p = reshape( U4(:,:,k-1) ,[],1);
    U5p = reshape( U5(:,:,k-1) ,[],1);
%     U6p = reshape( U6(:,:,k-1) ,[],1);
%     U7p = reshape( U7(:,:,k-1) ,[],1);
%     U8p = reshape( U8(:,:,k-1) ,[],1);
%     U9p = reshape( U9(:,:,k-1) ,[],1);
%     U10p = reshape( U10(:,:,k-1) ,[],1);
    


    % Get optimal sequence for player 1
    St1 = Sx{1,1}*x1k - Sy{1,1}*Yb + Sx{1,2}*x2k - Sy{1,2}*Yb + Su{1,2}*U2p;
    [U1o,J1o,exitflag] = quadprog(Rt{1},St1,M1,wr1);
    if exitflag<0
        error('Problems in the Optimization problem (player 1).');
    end
    
    % Get optimal sequence for player 2
    St2 = Sx{2,1}*x1k - Sy{2,1}*Yb + Sx{2,2}*x2k - Sy{2,2}*Yb + Su{2,1}*U1p ;
    [U2o,J2o,exitflag] = quadprog(Rt{2},St2,Mi,wr2);
    if exitflag<0
        error('Problems in the Optimization problem (player 2).');
    end
    
    % Get optimal sequence for player 3
    St3 = Sx{3,2}*x2k - Sy{3,2}*Yb + Sx{3,3}*x3k - Sy{3,3}*Yb + Su{3,2}*U2p ;
    [U3o,J3o,exitflag] = quadprog(Rt{3},St3,Mi,wr3);
    if exitflag<0
        error('Problems in the Optimization problem (player 3).');
    end
    
%     % Get optimal sequence for player 4
    St4 = Sx{4,3}*x3k - Sy{4,3}*Yb + Sx{4,4}*x4k - Sy{4,4}*Yb + Su{4,3}*U3p ;
    [U4o,J4o,exitflag] = quadprog(Rt{4},St4,Mi,wr4);
    if exitflag<0
        error('Problems in the Optimization problem (player 4).');
    end
    
    % Get optimal sequence for player 5
    St5 = Sx{5,4}*x4k - Sy{5,4}*Yb + Sx{5,5}*x5k - Sy{5,5}*Yb + Su{5,4}*U4p ;
    [U5o,J5o,exitflag] = quadprog(Rt{5},St5,Mi,wr5);
    if exitflag<0
        error('Problems in the Optimization problem (player 5).');
    end
    
%     % Get optimal sequence for player 6
%     St6 = Sx{6,5}*x5k - Sy{6,5}*Yb + Sx{6,6}*x6k - Sy{6,6}*Yb + Su{6,5}*U5p ;
%     [U6o,J6o,exitflag] = quadprog(Rt{6},St6,Mi,wr6);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 6).');
%     end
    
%     % Get optimal sequence for player 7
%     St7 = Sx{7,6}*x6k - Sy{7,6}*Yb + Sx{7,7}*x7k - Sy{7,7}*Yb + Su{7,6}*U6p ;
%     [U7o,J7o,exitflag] = quadprog(Rt{7},St7,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 7).');
%     end   
    
%     % Get optimal sequence for player 8
%     St8 = Sx{8,7}*x7k - Sy{8,7}*Yb + Sx{8,8}*x8k - Sy{8,8}*Yb + Su{8,7}*U7p ;
%     [U8o,J8o,exitflag] = quadprog(Rt{8},St8,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 8).');
%     end
%     
%     % Get optimal sequence for player 9
%     St9 = Sx{9,8}*x9k - Sy{9,8}*Yb + Sx{9,9}*x9k - Sy{9,9}*Yb + Su{9,8}*U8p ;
%     [U9o,J9o,exitflag] = quadprog(Rt{9},St9,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 9).');
%     end    
%   
%     % Get optimal sequence for player 10
%     St10 = Sx{10,9}*x10k - Sy{10,9}*Yb + Sx{10,10}*x10k - Sy{10,10}*Yb + Su{10,9}*U9p ;
%     [U10o,J10o,exitflag] = quadprog(Rt{10},St10,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 9).');
%     end
    
    
    
    
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
    U3pp = U3p;
    U4pp = U4p;
    U5pp = U5p;
%     U6pp = U6p;
%     U7pp = U7p;
%     U8pp = U8p;
%     U9pp = U9p;
%     U10pp = U10p;
    
    
    for pp = 1:np
        U1pp = w(1)*U1o + (1-w(1))*U1pp;
        U2pp = w(2)*U2o + (1-w(2))*U2pp;
        U3pp = w(3)*U3o + (1-w(3))*U3pp;
        U4pp = w(4)*U4o + (1-w(4))*U4pp;
        U5pp = w(5)*U5o + (1-w(5))*U5pp;
%         U6pp = w(6)*U6o + (1-w(6))*U6pp;
%         U7pp = w(7)*U7o + (1-w(7))*U7pp;
%         U8pp = w(8)*U8o + (1-w(8))*U3pp;
%         U9pp = w(8)*U9o + (1-w(9))*U9pp;
%         U10pp = w(10)*U10o + (1-w(10))*U10pp;
    end
    U1(:,:,k) = reshape( U1pp ,nu1,N);
    U2(:,:,k) = reshape( U2pp ,nu2,N);
    U3(:,:,k) = reshape( U3pp ,nu2,N);
    U4(:,:,k) = reshape( U4pp ,nu2,N);
    U5(:,:,k) = reshape( U5pp ,nu2,N);
%     U6(:,:,k) = reshape( U6pp ,nu2,N);
%     U7(:,:,k) = reshape( U7pp ,nu2,N);
%     U8(:,:,k) = reshape( U8pp ,nu2,N);
%     U9(:,:,k) = reshape( U9pp ,nu2,N);
%     U10(:,:,k) = reshape( U10pp ,nu2,N);
    
    % apply first value at each player
    u1k = U1(:,1,k);
    u2k = U2(:,1,k);
    u3k = U3(:,1,k);
    u4k = U4(:,1,k);
    u5k = U5(:,1,k);
%     u6k = U6(:,1,k);
%     u7k = U7(:,1,k);
%     u8k = U8(:,1,k);
%     u9k = U9(:,1,k);
%     u10k = U10(:,1,k);
    u(:,k) = [u1k+u(1,k-1)
              u2k+u(2,k-1) 
              u3k+u(3,k-1)
              u4k+u(4,k-1)
              u5k+u(5,k-1)
%               u6k+u(6,k-1)
%               u7k+u(7,k-1)
%               u8k+u(8,k-1)
%               u9k+u(9,k-1)
%               u10k+u(10,k-1)
                             ];
    W = (1/lambda - 1) * pk;
    % simulate system for distributed MPC
    %Xd(:,k+1) = A*Xd(:,k) + B*u(:,k) + Bd*dist(:,k); %simulate joint system
    %Xd(:,k+1) = A*Xd(:,k) + B*u(:,k);
    if k < 5
        Xd(:,k+1) = A*Xd(:,k) + B*u(:,k);
    else
        Xd(:,k+1) = A*Xd(:,k) + B(:,1)*u(1,k) + B(:,2)*u(2,k-1) + B(:,3)*u(3,k-2) + B(:,4)*u(4,k-3)+ B(:,5)*u(5,k-4); %delay
    end
    Y(:,k+1) = C*Xd(:,k+1) + noise(:,k);
    
%     %KALMAN FILTER
%     pk = A*pk*A' + W;
%     KalmanG = pk*C'*(C*pk*C'+V)^-1;
%     Xd(:,k+1) = Xd(:,k+1) + KalmanG*(Y(:,k+1)-C*Xd(:,k+1));
%     pk = (eye(numel(Xd(:,k+1)))-KalmanG*C) * pk;
end

%%
% figure(7101);
% 
% 
% grid on;
% hold on;
% plot(Xd(1,:),Xd(2,:),'-','Color',sstblue);
% plot(Xd(3,:), Xd(4,:), '-','Color',sstgray);
% % plot(Xd(5,:), Xd(6,:), '-','Color',sstgreen);
% % plot(Xd(7,:), Xd(8,:), '-','Color','red');
% % plot(Xd(9,:), Xd(10,:), '-','Color','magenta');
% 
% hold off;
% xlabel('$$x_1$$');
% ylabel('$$x_2$$');
% legend('car 1','car 2', 'car 3', 'car 4', 'car 5', 'car 6', 'car 7', 'car 8', 'car 9', 'car 10');
% title('Phase plot');

figure();
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(TX,Y(1,:),'-','Color',sstblue);
% plot(TX,Xd(2,:),'-','Color',sstgray);
plot(TX,Y(2,:),'-','Color',sstgreen);
% plot(TX,Xd(4,:),'-','Color','magenta');
plot(TX,Y(3,:),'-','Color','red');
% % plot(TX,Xd(6,:),'-','Color','cyan');
plot(TX,Y(4,:),'-','Color','magenta');
% % plot(TX,Xd(8,:),'-','Color','yellow');
plot(TX,Y(5,:),'-','Color','green');
% % plot(TX,Xd(10,:),'-','Color','blue');
% plot(TX,Xd(11,:),'--','Color',sstblue);
% % plot(TX,Xd(12,:),'--','Color',sstgray);
% plot(TX,Xd(13,:),'--','Color',sstgreen);
% plot(TX,Xd(14,:),'--','Color','magenta');

hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r$$','car1 $$x_1$$ ','car1 $$x_2$$','car2$$x_1$$','car2 $$x_2$$', 'car3 $$x_1$$ ','car3 $$x_2$$','car4$$x_1$$','car4 $$x_2$$', 'car5 $$x_1$$ ','car5 $$x_2$$','car 6$$x_1$$','car6 $$x_2$$', 'car7 $$x_1$$ ','car7 $$x_2$$','car8$$x_1$$','car8 $$x_2$$', 'car9 $$x_1$$ ','car9 $$x_2$$','car10$$x_1$$','car10 $$x_2$$','Location','SouthEast');
title('State evolution');

% figure(7103);
% plot(TU,u(1,:),'d-','Color',sstblue);
% grid on;
% hold on;
% plot(TU,u(2,:),'d-','Color',sstgray);
% plot(TU,u(3,:),'d-','Color',sstgreen);
% plot(TU,u(4,:),'-','Color','magenta');
% % plot(TU,u(5,:),'-','Color','yellow');
% % plot(TU,u(6,:),'-','Color','green');
% % plot(TU,u(7,:),'-','Color','blue');
% 
% hold off;
% 
% 
% xlabel('$$t_k$$');
% ylabel('$$u(t_k)$$');
% legend('$$u_1$$','$$u_2$$','$$u_3$$','$$u_4$$','$$u_5$$','$$u_6$$','$$u_7$$','$$u_8$$','$$u_9$$','$$u_{10}$$');
% title('Input');










