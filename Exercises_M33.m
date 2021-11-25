%% Exercise book solutions (CPCS)
%  Module 3: Nonlinear model predictive control and stability analysis
%
% Bruno Guerreiro (bj.guerreiro@fct.unl.pt)

clear all;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstgray         = [70,70,70]/255;
sstlightgray    = [200,200,200]/255;

% Check if Multi-parametric-toolbox is installed (from https://www.mpt3.org/Main/Installation)
try
    Polyhedron(0,0);
    disp('Great! Multi-parametric-toolbox seems to be installed.');
catch
    error('Please install Multi-parametric-toolbox (from https://www.mpt3.org/Main/Installation).');
end

%% Exercise 3.1: simple quasi-linear optimization

% specify system parameters
a = 0.99;
b = 0.1;
c = 0.1;

x_0 = 0;

% performance index parameters
r1 = 0.001;
r2 = 0.1;
r3 = 1;
r4 = 10;
N = 99;
r_N = 10;

% compute state, costate and control evolution
[x1,u1,lambda1,K] = OptCtrl_1(a,b,c,r1,N,x_0,r_N);
[x2,u2,lambda2,K] = OptCtrl_1(a,b,c,r2,N,x_0,r_N);
[x3,u3,lambda3,K] = OptCtrl_1(a,b,c,r3,N,x_0,r_N);
[x4,u4,lambda4,K] = OptCtrl_1(a,b,c,r4,N,x_0,r_N);

figure(5101);
subplot(311);
stairs(K,x1,'-','Color',sstdarkblue);
hold on;
stairs(K,x2,'--','Color',sstblue);
stairs(K,x3,'-.','Color',sstlightblue);
stairs(K,x4,'--','Color',sstgreen);
stairs(K,r_N*ones(size(x1)),'-.','Color',sstlightgreen);
hold off;
grid on;
xlabel('$$k$$');
ylabel('$$x(k)$$');
legend(['$$r = ' num2str(r1) '$$'],['$$r = ' num2str(r2) '$$'],['$$r = ' num2str(r3) '$$'],['$$r = ' num2str(r4) '$$'],'reference','Location','west');
subplot(312);
stairs(K,lambda1,'-','Color',sstdarkblue);
hold on;
stairs(K,lambda2,'--','Color',sstblue);
stairs(K,lambda3,'-.','Color',sstlightblue);
stairs(K,lambda4,'--','Color',sstgreen);
hold off;
grid on;
xlabel('$$k$$');
ylabel('$$\lambda(k)$$');
subplot(313);
stairs(K,u1,'-','Color',sstdarkblue);
hold on;
stairs(K,u2,'--','Color',sstblue);
stairs(K,u3,'-.','Color',sstlightblue);
stairs(K,u4,'--','Color',sstgreen);
hold off;
grid on;
xlabel('$$k$$');
ylabel('$$u(k)$$');

%% Exercise 3.2: NMPC for car on horizontal plane

clear X U Xopt Uopt TX TU;
example_name = 'simpleNLcar_MPC_simul';

% Define the model parameters:
%MParam.beta = 0.5;
%MParam.m = 1;
MParam.Ts = 0.1;

% define the cost parameters
N = 5;
Q = diag([8,0.2]);
P = 10*Q;
R = 0.3;
x0 = [0.7;-0.7];
u0 = 0;
nx = size(x0,1);
nu = size(u0,2);

% compute the cost batch matrices
[~,~,Qb,Rb] = GetBatchXMatrices(zeros(nx),zeros(nx,nu),[],N,P,Q,R);

%simulate controlled system:
nk = 60;
TU = 1:nk;
TX = 1:nk+1;
X = kron(x0,ones(1,nk+1));
U = kron(u0,ones(1,nk));
Xopt(:,:,1) = X(:,1:N+1);
Uopt(:,:,1) = U(:,1:N);
for k = 1:nk
    TUopt(k,:) = k:k+N-1;
    TXopt(k,:) = k:k+N;
    
    % inicialize optimization problem based on previous solution:
    xk = X(:,k);
    Xopt(:,1,k) = xk;
    X0 = reshape( Xopt(:,:,k) , [],1);
    U0 = reshape( Uopt(:,:,k) , [],1);
    XU0 = [X0;U0];
    
    % define model constraint and cost function handlers:
    fun = @(XU)	cost(XU,Qb,Rb);
    nonlcon = @(XU)	statecon(XU,xk,N,nx,nu,MParam);
    [XUopt,Jo,exitflag,output] = fmincon(fun,XU0,[],[],[],[],[],[],nonlcon);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end
    Xopt(:,:,k) = reshape( XUopt(1 : nx*(N+1))               , nx,[]);
    Uopt(:,:,k) = reshape( XUopt(nx*(N+1)+1:nx*(N+1)+nu*N)   , nu,[]);
    
    % get MPC policy and simulate the nonlinear system:
    U(:,k) = Uopt(:,1,k);
    X(:,k+1) = f_pend(X(:,k),U(:,k),MParam);
    
    % next initialization
    Xopt(:,1:N,k+1) = Xopt(:,2:N+1,k);
    Xopt(:,N+1,k+1) = Xopt(:,N+1,k);
    Uopt(:,1:N-1,k+1) = Uopt(:,2:N,k);
    Uopt(:,N,k+1)   = Uopt(:,N,k);
    
end

figure(5201);
plot(X(1,:),X(2,:),'s-','Color',sstblue);
grid on;
hold on;
for k = 1:nk
    plot(Xopt(1,:,k),Xopt(2,:,k),'.--','Color',sstlightblue);
end
plot(X(1,:),X(2,:),'s-','Color',sstblue);
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('const.','pred.');
title('Phase plot');

figure(5202);
plot(TX,X(1,:),'s-','Color',sstblue);
grid on;
hold on;
plot(TX,X(2,:),'d-','Color',sstgreen);
for k = 1:nk
    plot(TXopt(k,:),Xopt(1,:,k),'.--','Color',sstlightblue);
    plot(TXopt(k,:),Xopt(2,:,k),'.-.','Color',sstlightgreen);
end
plot(TX,X(1,:),'s-','Color',sstblue);
plot(TX,X(2,:),'d-','Color',sstgreen);
hold off;
xlabel('$$t_k$$ [s]');
ylabel('$$x(t_k)$$');
title('State evolution');
legend('$$x_1$$ const.','$$x_2$$ const.','$$x_1$$ pred.','$$x_2$$ pred.','Interpreter','latex');

figure(5203);
plot(TU,U,'s-','Color',sstblue);
grid on;
hold on;
for k = 1:nk
    plot(TUopt(k,:),Uopt(:,:,k),'.--','Color',sstlightblue);
end
plot(TU,U,'s-','Color',sstblue);
hold off;
xlabel('$$t_k$$ [s]');
ylabel('$$u(t_k)$$');
legend('const.','pred.');
title('Input');

%% Exercise 3.3: no solution provided (mass-spring-dumper NMPC)

%% Exercise 3.4: one-step controllable set

A = [   0.5  0
        1   -0.5 ];
x_max = [10;10];
H = [-eye(2);eye(2)];
h = 10*ones(4,1);

X = Polyhedron(H,h); % or X = Polyhedron('lb',-x_max,'ub',x_max);
X.computeVRep();
% setX_.computeVRep();

sys = LTISystem('A', A);
sys.x.min = -x_max;
sys.x.max = x_max;
PreX = sys.reachableSet('Direction','backward');
K1 = PreX.intersect(X).minHRep();

X.H,
PreX.H,
K1.H,

figure(6101);
plot(X,'color',sstblue,PreX,'color',sstlightblue,K1,'color',sstgreen);
xlabel('$$x_1$$','Interpreter','Latex');
ylabel('$$x_2$$','Interpreter','Latex');
legend('$$\mathcal{X}$$','Pre($$\mathcal{X}$$)','$$\mathcal{K}_1(\mathcal{X})$$','Interpreter','Latex','Location','Southeast');

%% Exercise 3.5: maximal controllable set

% define the system
A = [   1.5 0
        1   -1.5 ];
B = [1;0];
sys = LTISystem('A',A,'B',B);

% define the model constraint sets
x_max = [10;10];
u_max = 5;
sys.x.min = -x_max;
sys.x.max = x_max;
sys.u.min = -u_max;
sys.u.max = u_max;

% % Sets X and U defined within system, but can be obtained independently:
% X = Polyhedron('lb',-x_max,'ub',x_max);
% U = Polyhedron('lb',-u_max,'ub',u_max);
% X = sys.x.boundsToPolyhedron();
% U = sys.u.boundsToPolyhedron();

% define target set:
S = Polyhedron('lb',-ones(2,1),'ub',ones(2,1));

% compute some n-step controllable sets and the maximal
K1 = GetKn(S,sys,1);
K2 = GetKn(S,sys,2);
K3 = GetKn(S,sys,3);
[Kmax,isf,Nb] = GetKn(S,sys,1000);

% plot resulting sets
figure(6202);
plot(Kmax,'color',sstblue,K3,'color',sstgreen,K2,'color',sstlightblue,K1,'color',sstlightgreen,S,'color',sstgray);
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend( '$$\mathcal{K}_{max}(\mathcal{S})$$',...
        '$$\mathcal{K}_3(\mathcal{S})$$',...
        '$$\mathcal{K}_2(\mathcal{S})$$',...
        '$$\mathcal{K}_1(\mathcal{S})$$',...
        '$$\mathcal{S}$$',...
        'Location','Southeast');

%% Exercise 3.6: no solution provided (maximal controllable set)

%% Exercise 3.7: MPC stability analysis

% define system
A = [1.4 1;0 0.3];
B = [0;0.5];
sys = LTISystem('A',A,'B',B);

% define the model constraint sets
x_max = [12;12];
u_max = 3;
sys.x.min = -x_max;
sys.x.max = x_max;
sys.u.min = -u_max;
sys.u.max = u_max;
X = Polyhedron('lb',-x_max,'ub',x_max);
U = Polyhedron('lb',-u_max,'ub',u_max);

% Case (a): define horizon and final constraint set X_f = {x = 0}:
N = 2;
Xf = Polyhedron('Ae', eye(2), 'be', zeros(2,1));

% get the feasibility set X0:
[X0a,isfinite,Nb,Xi] = GetXi(Xf,sys,N);

% plot resulting sets
figure(6301);
plot(X0a,'color',sstblue,Xi(end-1),'color',sstgreen,Xf,'color',sstlightblue);
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend( '$$\mathcal{X}_0$$',...
        '$$\mathcal{X}_1$$',...
        '$$\mathcal{X}_f$$');
title('Feasibility set, case a)');
    
% Case (b): define horizon and final constraint set X_f = {x = 0}:
N = 42;
Q = eye(2);
R = 0.1;
sys.x.penalty = QuadFunction(Q);
sys.u.penalty = QuadFunction(R);
Xf = sys.LQRSet();

% get the feasibility set X0:
[X0b,isfinite,Nb,Xi] = GetXi(Xf,sys,N);

% plot resulting sets
figure(6302);
plot(X0b,'color',sstblue,Xi(end-1),'color',sstgreen,Xf,'color',sstlightblue);
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend( '$$\mathcal{X}_0$$',...
        '$$\mathcal{X}_{N-1}$$',...
        '$$\mathcal{X}_f$$');
title('Feasibility set, case b)');

% Case (c): define horizon and final constraint set X_f = R^2:
xf_max = 1e3*ones(2,1);
Xf = Polyhedron('lb',-xf_max,'ub',xf_max);
[~,isf,Nb] = GetKn(Xf,sys,1000);
N = 42; % N = Nb + 1;

% get the feasibility set X0:
[X0c,isfi,Nbi,Xi] = GetXi(Xf,sys,N);

% plot resulting sets
figure(6303);
plot(Xi(end-1),'color',sstgreen,Xi(end-2),'color',sstgray,Xi(end-3),'color',sstlightblue,...
     Xi(end-4),'color',sstlightgreen,Xi(end-5),'color',sstlightgray,X0c,'color',sstblue);
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend( '$$\mathcal{X}_{N-1}$$',...
        '$$\mathcal{X}_{N-2}$$',...
        '$$\mathcal{X}_{N-3}$$',...
        '$$\mathcal{X}_{N-4}$$',...
        '$$\mathcal{X}_{N-5}$$',...
        '$$\mathcal{X}_0$$');
title('Feasibility set');

%% Exercise 3.8: no solution provided (MPC stability analysis)


%% nested functions %%

% define cost function
function J = cost(XU,Qb,Rb)
    J = XU'*blkdiag(Qb,Rb)*XU;
end

% constraints for simple car:
function [c,ceq] = statecon(XU,x0,N,nx,nu,MP)
    X = reshape( XU(1 : nx*(N+1))               , nx,[]);
    U = reshape( XU(nx*(N+1)+1:nx*(N+1)+nu*N)   , nu,[]);
    Ceq = zeros(nx,N);
    Ceq(:,1) = X(:,1) - x0;
    for k = 1:N
        Ceq(:,k+1) = f_pend(X(:,k),U(:,k),MP) - X(:,k+1);
    end
    ceq = reshape(Ceq,[],1); c = [];
end

function xp = f_car(x,u,P)
    xdot = [x(2);1/P.m*(u(1)-P.beta*x(2)^2)];
    xp = x + P.Ts*xdot;
end

function xp = f_pend(x,u,P)
    xdot = [x(2); u(1) - 0.2*x(2)^2 - 10*sin(x(1))];
    xp = x + P.Ts*xdot;
end
