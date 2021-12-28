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
%% setup & model
m = 1500;
Cd = 0.3;
area = 1.8;
rho = 1.225;
g = 9.8065;
theta = 0.05;
fmax = 2800;
fmin = -2200;
deltaFmax = 200;%300;
maxTheta = 0.1;
dDist = 10;
distMin = 3;
vRe = 25;
ve = 25;
we = 0;
Ts = 0.1;
d = (rho*area*Cd*ve)/m;

% Ac = [0  -1  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
%      0  -d  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
%      0  1   0   -1  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
%      0  0   0   -d  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
%      0  0   0   1   0   -1  0   0   0   0   0   0   0   0   0   0   0   0   0   0
%      0  0   0   0   0   -d  0   0   0   0   0   0   0   0   0   0   0   0   0   0
%      0  0   0   0   0   1   0   -1  0   0   0   0   0   0   0   0   0   0   0   0
%      0  0   0   0   0   0   0   -d  0   0   0   0   0   0   0   0   0   0   0   0
%      0  0   0   0   0   0   0   1   0   -1  0   0   0   0   0   0   0   0   0   0
%      0  0   0   0   0   0   0   0   0   -d  0   0   0   0   0   0   0   0   0   0
%      0  0   0   0   0   0   0   0   0   1   0   -1  0   0   0   0   0   0   0   0
%      0  0   0   0   0   0   0   0   0   0   0   -d  0   0   0   0   0   0   0   0
%      0  0   0   0   0   0   0   0   0   0   0   1   0   -1  0   0   0   0   0   0
%      0  0   0   0   0   0   0   0   0   0   0   0   0   -d  0   0   0   0   0   0
%      0  0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1  0   0   0   0
%      0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   -d  0   0   0   0
%      0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1  0   0
%      0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -d  0   0
%      0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1
%      0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -d];
 Ac = [0  -1  0   0        
       0  -d  0   0           
       0  1   0   -1         
       0  0   0   -d      ];    
% Bc= [  0   0    0   0   0   0   0   0   0   0
%        1/m 0    0   0   0   0   0   0   0   0
%         0   0   0   0   0   0   0   0   0   0
%         0  1/m  0   0   0   0   0   0   0   0
%         0   0   0   0   0   0   0   0   0   0
%         0   0   1/m 0   0   0   0   0   0   0
%         0   0   0   0   0   0   0   0   0   0
%         0   0   0   1/m 0   0   0   0   0   0
%         0   0   0   0   0   0   0   0   0   0
%         0   0   0   0   1/m 0   0   0   0   0  
%         0   0   0   0   0   0   0   0   0   0
%         0   0   0   0   0   1/m 0   0   0   0
%         0   0   0   0   0   0   0   0   0   0
%         0   0   0   0   0   0   1/m 0   0   0
%         0   0   0   0   0   0   0   0   0   0
%         0   0   0   0   0   0   0   1/m 0   0
%         0   0   0   0   0   0   0   0   0   0
%         0   0   0   0   0   0   0   0   1/m 0
%         0   0   0   0   0   0   0   0   0   0
%         0   0   0   0   0   0   0   0   0   1/m ];
Bc= [  0   0       
       1/m 0          
        0   0        
        0  1/m       

 ];
% Bijc = {Bc(1:2,1),Bc(1:2,2),Bc(1:2,3),Bc(1:2,4),Bc(1:2,5),Bc(1:2,6),Bc(1:2,7),Bc(1:2,8),Bc(1:2,9),Bc(1:2,10)
%         Bc(2:4,1),Bc(2:4,2),Bc(2:4,3),Bc(2:4,4),Bc(2:4,5),Bc(2:4,6),Bc(2:4,7),Bc(2:4,8),Bc(2:4,8),Bc(2:4,10)
%         Bc(4:6,1),Bc(4:6,2),Bc(4:6,3),Bc(4:6,4),Bc(4:6,5),Bc(4:6,6),Bc(4:6,7),Bc(4:6,8),Bc(4:6,9),Bc(4:6,10)
%         Bc(6:8,1),Bc(6:8,2),Bc(6:8,3),Bc(6:8,4),Bc(6:8,5),Bc(6:8,6),Bc(6:8,7),Bc(6:8,8),Bc(6:8,9),Bc(6:8,10)
%         Bc(8:10,1),Bc(8:10,2),Bc(8:10,3),Bc(8:10,4),Bc(8:10,5),Bc(8:10,6),Bc(8:10,7),Bc(8:10,8),Bc(8:10,9),Bc(8:10,10)
%         Bc(10:12,1),Bc(10:12,2),Bc(10:12,3),Bc(10:12,4),Bc(10:12,5),Bc(10:12,6),Bc(10:12,7),Bc(10:12,8),Bc(10:12,9),Bc(10:12,10)
%         Bc(12:14,1),Bc(12:14,2),Bc(12:14,3),Bc(12:14,4),Bc(12:14,5),Bc(12:14,6),Bc(12:14,7),Bc(12:14,8),Bc(12:14,9),Bc(12:14,10)
%         Bc(14:16,1),Bc(14:16,2),Bc(14:16,3),Bc(14:16,4),Bc(14:16,5),Bc(14:16,6),Bc(14:16,7),Bc(14:16,8),Bc(14:16,9),Bc(14:16,10)
%         Bc(16:18,1),Bc(16:18,2),Bc(16:18,3),Bc(16:18,4),Bc(16:18,5),Bc(16:18,6),Bc(16:18,7),Bc(16:18,8),Bc(16:18,9),Bc(16:18,10)
%         Bc(18:20,1),Bc(18:20,2),Bc(18:20,3),Bc(18:20,4),Bc(18:20,5),Bc(18:20,6),Bc(18:20,7),Bc(18:20,8),Bc(18:20,9),Bc(18:20,10)};

Bijc = {Bc(1:2,1),Bc(1:2,2)
        Bc(2:4,1),Bc(2:4,2)
        
       };
% C=[ 1  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
%     0  0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
%     0  0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
%     0  0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0
%     0  0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0
%     0  0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0
%     0  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0
%     0  0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0
%     0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0
%     0  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0];
C=[ 1  0   0   0       
    0  0   1   0     ];
% Ci={[1 0], [0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0],[0 1 0]};

 Ci={[1 0], [0 1 0]};
 A1c = [0 -1;0 -d];
 A2c = [-d 0 0;1 0 -1;0 0 -d];
 
% A = eye(20) + Ac*Ts;
% Ai={eye(2) + A1c*Ts,eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts,...
%      eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts, eye(3) + A2c*Ts,...
%      eye(3) + A2c*Ts,eye(3) + A2c*Ts};
A=eye(4) + Ac*Ts;
Ai={eye(2) + A1c*Ts, eye(3) + A2c*Ts};
B = Bc*Ts;
for i=1:2
    for j=1:2
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
xd0 = ones(4,1);


%A_ext = [A zeros(nx,ny);C*A eye(ny)];
Ai_ext{1} = [Ai{1} zeros(nx1,ny1);Ci{1}*Ai{1} eye(ny1)];
for i=2:2
   Ai_ext{i} = [Ai{i} zeros(nx2,ny2);Ci{i}*Ai{i} eye(ny2)];
end

for i=1:2
   Bij_ext{1,i} = [Bij{1,i};Ci{1}*Bij{1,i}]; 
end
 
for i=2:2
   for j = 1:2
        Bij_ext{i,j} = [Bij{i,j};Ci{i}*Bij{i,1}]; 
   end
end
 
Ci_ext{1} = [zeros(ny1,nx1),eye(ny1)];

for i=2:2
   Ci_ext{i} = [zeros(ny2,nx2),eye(ny2)]; 
end
 
 %% Cost
 N = 60;
Pi = 150;
Qi = 5000;
Ri = 10;
alphai = [1 1 1 1 1 1 1 1 1 1 ];
% distributed steps parameters
w(1) = 0.5;
w(2) = 0.5;
% for i=2:3
%    w(i) = (1-w(1))/9; 
% end
np=2;

%% Batch matrices
[Fb,Gb,Qb,Rb] = GetBatchYNdistMatrices(Ai_ext,Bij_ext,Ci_ext,N,Pi,Qi,Ri,alphai(1));
[Rt,Sx,Sy,Su,K,Ky,L] = GetYMPC(Fb,Gb,Qb,Rb);

% %% Constraints 
% 
% %player 1
% 
% u_max1 = fmax -(0.5*rho*area*Cd*ve^2);
% u_min1 = fmin - (0.5*rho*area*Cd*ve^2);
% U_max1 = kron(u_max1,ones(N*nu1,1));
% U_min1 = kron(u_min1,ones(N*nu1,1));
% M31 = kron(tril(ones(N)), eye(nu1));
% M41 = kron(ones(N,1), eye(nu1));
% Mu1 = [-M31;M31];
% 
% du_max1 = deltaFmax;
% du_min1 = -deltaFmax;
% DU_max1 = kron(du_max1,ones(N*nu1,1));
% DU_min1 = kron(du_min1,ones(N*nu1,1));
% Mdu1 = [-eye(N*nu1);eye(N*nu1)];
% wdu1 = [-DU_min1;DU_max1];
% 
% pr_min = -7 ;
% pr_max = 93;
% Y_min1 = kron(pr_min,ones(nu1*(N+1),1));
% Y_max1 = kron(pr_max, ones(nu1*(N+1),1));
% My1 = [-Gb11; Gb11];
% 
% M1 = [Mdu1; Mu1];
% 
% 
% 
% 
% 
% %player i
% 
% u_max2 = fmax -(0.5*rho*area*Cd*ve^2);
% u_min2 = fmin - (0.5*rho*area*Cd*ve^2);
% U_max2 = kron(u_max2,ones(N*nu2,1));
% U_min2 = kron(u_min2,ones(N*nu2,1));
% M32 = kron(tril(ones(N)), eye(nu2));
% M42 = kron(ones(N,1), eye(nu2));
% Mu2 = [-M32;M32];
% 
% du_max2 = deltaFmax;
% du_min2 = -deltaFmax;
% DU_max2 = kron(du_max2,ones(N*nu2,1));
% DU_min2 = kron(du_min2,ones(N*nu2,1));
% Mdu2 = [-eye(N*nu2);eye(N*nu2)];
% wdu2 = [-DU_min2;DU_max2];
% 
% pr_min = -7 ;
% pr_max = 93;
% Y_min2 = kron(pr_min,ones(nu2*(N+1),1));
% Y_max2 = kron(pr_max, ones(nu2*(N+1),1));
% My2 = [-Gb22; Gb22];
% 
% M2 = [Mdu2; Mu2];














