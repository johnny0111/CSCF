%% Exercise book solutions (CPCS)
%  Module 5: Hybrid MPC
%
% Bruno Guerreiro (bj.guerreiro@fct.unl.pt)

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

% Check if Multi-parametric-toolbox is installed (from https://www.mpt3.org/Main/Installation)
try
    Polyhedron(0,0);
    disp('Great! Multi-parametric-toolbox seems to be installed.');
catch
    error('Please install Multi-parametric-toolbox (from https://www.mpt3.org/Main/Installation).');
end

%% Exercise 5.1: boucing ball hybrid system

%bouncing ball
S.Ts=0.01;
S.g=9.8;
S.alpha=0.3;
S.h_min = 0;
S.h_max = 6;
S.v_max = 100;
model = MLDSystem('bouncing_ball',S);

Tsim = 7;
Nsim=Tsim/S.Ts;
U=zeros(0,Nsim);
x0=[5 0]';
model.initialize(x0);
data = model.simulate(U);
X = data.X;
TX = 0:S.Ts:S.Ts*Nsim;
trigger = (X(1,:)<=0)+1;

figure(8100);
A1 = [-eye(2);eye(2)];
b1 = [eps;100;10;100];
R1 = Polyhedron(A1,b1);
A2 = [-eye(2);eye(2)];
b2 = [10;100;0;100];
R2 = Polyhedron(A2,b2);
plot(R1,'color',sstlightblue,'FaceAlpha',0.1,R2,'color',sstgreen,'FaceAlpha',0.1);
xlabel('$$h(k)$$');
ylabel('$$v(k)$$');
title('PWA regions');
text(5,-80,'$$\mathcal{R}_1: h(k) \geq \epsilon$$','HorizontalAlignment','center','FontSize',20)
text(-5,-80,'$$\mathcal{R}_2: h(k) \leq 0$$','HorizontalAlignment','center','FontSize',20)

figure(8101);
plot(X(1,:),X(2,:),'.-','Color',sstblue);
grid on;
xlabel('$$h(k)$$');
ylabel('$$v(k)$$');
title('Phase plot');

figure(8102);
stairs(TX,X(1,:),'.-','Color',sstblue);
grid on;
hold on;
stairs(TX,X(2,:),'.--','Color',sstgreen);
stairs(TX,trigger,'.-.','Color',sstlightblue);
hold off;
xlabel('$$t_k$$ [s]');
title('State evolution');
legend('$$h(k)$$','$$v(k)$$','$$i(k)$$');

%% Exercise 4.2: no solution provided (2-D Bouncing ball)

%% Exercise 4.3: basic distributed tracking MPC

% simple pendulum (inverted or not)
% nonlinear equation: I dot(om) = tau + m g l sin(th) - beta om
S.m = 1;
S.l = 0.8;
S.g = 9.81;
S.beta = 1.5;
S.I = S.m*S.l^2;
S.Ts = 0.05;
S.tau_min = 0.1;
S.tau_max = 11;
S.th_max = 2*pi;
S.om_max = 20;

% nonlinear function piecewise approximation: 
% sin(th) = { alpha th IF |th| <= pi/2 ; -alpha th - gamma IF th <= -pi/2 ; -alpha th + gamma IF th >= pi/2  }
S.alpha = 24/pi^3;
S.gamma = 24/pi^2;

% testing
TH = -2*pi:0.1:2*pi;
f1 = sin(TH);
f2 = approx_sin(TH,S);
figure(8300);
plot(TH,f1,'g.-.');
grid on;
hold on;
plot(TH,f2,'b.--');
hold off;
xlabel('$$\theta$$');
legend('$$\sin(\theta)$$','PWA approx. $$\sin(\theta)$$');
title('Sin approximation');

% create hybrid pendulum model
Ac = [0, 1; S.g/S.l*S.alpha, -S.beta/S.I ];
Bc = [0,0;S.g/S.l,1/S.I];
Cc = [1,0];
Dc = zeros(1,2);
nx = 2;
% exact discretization
lin_pend = ss(Ac,Bc,Cc,Dc);
disc_pend = c2d(lin_pend,S.Ts);
[Ad2,Bd2,Cd2,Dd2] = ssdata(disc_pend);
% forward Euler discretization
Ad = eye(nx) + S.Ts*Ac;
Bd = S.Ts*Bc;
S.a11 = Ad(1,1);
S.a12 = Ad(1,2);
S.a21 = Ad(2,1);
S.a22 = Ad(2,2);
S.b11 = Bd(1,1);
S.b12 = Bd(1,2);
S.b21 = Bd(2,1);
S.b22 = Bd(2,2);
model = MLDSystem('hyb_pendulum',S);

% Open loop simulation of hybrid and nonlinear models
N=500;
tX = 0:S.Ts:S.Ts*N;
tRef = 0:S.Ts:S.Ts*(N-1);
U=[zeros(1,round(N*3/5)),2*ones(1,round(N*2/5))];
x0=[0.01;0];
model.initialize(x0);
data = model.simulate(U);
X = data.X;
Xnl = simul_NL_pend(U,x0,S);

figure(8301);
plot(X(1,:),X(2,:),'.-','Color',sstblue);
grid on;
hold on;
plot(Xnl(1,:),Xnl(2,:),'.--','Color',sstgreen);
hold off;
xlabel('$$\theta$$');
ylabel('$$\dot{\theta}$$');
title('Phase plot');
legend('hybrid','nonlinear');

figure(8302);
plot(tRef,U,'k-');
grid on;
hold on;
plot(tX,X(1,:),'.--','Color',sstblue);
plot(tX,X(2,:),'.-.','Color',sstlightblue);
plot(tX,Xnl(1,:),'.--','Color',sstgreen);
plot(tX,Xnl(2,:),'.-.','Color',sstlightgreen);
hold off;
xlabel('t [s]');
title('State evolution');
legend('$$\tau$$','$$\theta_{hyb}$$','$$\dot{\theta}_{hyb}$$','$$\theta_{nl}$$','$$\dot{\theta}_{nl}$$');

% design hybrid mpc for above model
N = 5;
pred_model = model;
pred_model.u.min = -10;
pred_model.u.max = 10;
pred_model.x.with('reference');
pred_model.x.reference = 'free';
% pred_model.x.penalty = QuadFunction( diag([1, 0]) );
% pred_model.u.penalty = QuadFunction( 0.01 );
pred_model.x.penalty = OneNormFunction( diag([1, 0]) );
pred_model.u.penalty = OneNormFunction( 0.01 );
ctrl = MPCController(pred_model,N);

% simulate closed loop
x0 = [pi; 0];
Nsim = 100;
tX = 0:S.Ts:S.Ts*Nsim;
tRef = 0:S.Ts:S.Ts*(Nsim-1);
xref = [pi;0]*(tRef>=S.Ts*Nsim/2);
% loop = ClosedLoop(ctrl, model);
% dataCL = loop.simulate(x0, Nsim, 'x.reference', xref);
% Xcl = dataCL.X;
nx = length(x0);
Xcl = zeros(nx,Nsim+1);
Xcl(:,1) = x0;
Ucl = zeros(1,Nsim);
for k = 1:Nsim
   Ucl(:,k) = ctrl.evaluate(Xcl(:,k), 'x.reference', xref(:,k));
   Xcl(:,k+1) = f_pend(Xcl(:,k),Ucl(:,k),S);
   fprintf('HMPC closed-loop k = %d/%d\n',k,Nsim);
end

figure(8311);
plot(xref(1,:),xref(2,:),'k-');
grid on;
hold on;
plot(Xcl(1,:),Xcl(2,:),'.-','Color',sstblue);
hold off;
xlabel('$$\theta$$');
ylabel('$$\dot{\theta}$$');
title('Phase plot');
legend('ref','system');

figure(8312);
stairs(tRef,xref(1,:),'k-');
grid on;
hold on;
stairs(tX,Xcl(1,:),'.--','Color',sstblue);
stairs(tX,Xcl(2,:),'.-.','Color',sstlightblue);
hold off;
xlabel('t [s]');
title('State evolution');
legend('$$\tau$$','$$\theta_{hyb}$$','$$\dot{\theta}_{hyb}$$');

figure(8313);
stairs(tRef,Ucl(1,:),'.-','Color',sstblue);
grid on;
hold on;
plot(tRef,pred_model.u.min*ones(size(tRef)),'--','Color',sstgreen);
plot(tRef,pred_model.u.max*ones(size(tRef)),'-.','Color',sstlightgreen);
hold off;
xlabel('t [s]');
ylabel('$$\tau(t_k)$$');
title('Input');
legend('$$\bar{y}$$','$$\tau_{min}$$','$$\tau_{max}$$');

%% Exercise 5.4: no solution provided (2x2 PWA)

%% %%%%%%%%%%%%%%%%%% auxiliary functions %%%%%%%%%%%%%%%

function out = approx_sin(th,S)

    out = zeros(size(th));
    out(abs(th) <= pi/2) =  S.alpha*th(abs(th) <= pi/2);
    out(th < -pi/2)      = -S.alpha*th(th < -pi/2)   - S.gamma;
    out(th >  pi/2)      = -S.alpha*th(th >  pi/2)   + S.gamma;
    out(th < -pi*3/2)    =  S.alpha*th(th < -pi*3/2) + 2*S.gamma;
    out(th >  pi*3/2)    =  S.alpha*th(th >  pi*3/2) - 2*S.gamma;
    
end

function X = simul_NL_pend(U,x0,P)
    N = size(U,2);
    nx = length(x0);
    X = zeros(nx,N+1);
    X(:,1) = x0;
    for k = 1:N
        X(:,k+1) = f_pend(X(:,k),U(:,k),P);
    end
end

function xp = f_pend(x,u,P)
    % nonlinear equation: dot(th) = om; I dot(om) = u + m g l sin(th) - beta om
    xdot = [x(2);1/P.I*(u(1) + P.m*P.l*P.g*sin(x(1))) - P.beta*x(2)];
    xp = x + P.Ts*xdot;
end