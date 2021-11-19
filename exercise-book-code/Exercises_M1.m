%% Exercise book solutions (CPCS)
%  Module 1: Modeling, optimization and optimal control
%
% Bruno Guerreiro (bj.guerreiro@fct.unl.pt)

sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstgreen        = [43,191,92]/255;

%% Exercise 1.1

A = [0 1;-1 -4];
B = [0;0.5];
C = [1 0];
D = 0;

% by the definition:
s = tf('s');
Hc1 = C*(s*eye(2)-A)^-1*B+D,

% using matlab:
sys = ss(A,B,C,D);
Hc2 = tf(sys),

%% Exercise 1.2

A = 2;
B = 1;
C = 0.5;
D = 0;

% by the definition:
z = tf('z',1);
Hd1 = C*(z*eye(1)-A)^-1*B+D,

% using matlab:
sys = ss(A,B,C,D,1);
Hd2 = tf(sys),

%% Exercise 1.3

beta = 0.5;
m = 1;

Ac = [0 1;0 -beta/m];
Bc = [0;1/m];
Cc = eye(2);
Dc = zeros(2,1);
lbd_c = eig(Ac),

% simulate systems;
t = (0:0.01:10)';
% u = t>=1;
u = t*0;
x0 = [1;1];
sys_c = ss(Ac,Bc,Cc,Dc);
y_c = lsim(sys_c,u,t,x0);

% discretization:
T = 0.4;
Ad1 = expm(Ac*T);
Bd1 = T*Bc; % approximation (actual Bd1 might be time-varying)
Ad2 = eye(2) + T*Ac;
Bd2 = T*Bc;
Cd = Cc;
Dd = 0;
lbd_d = eig(Ad1),

td = (0:T:10)';
% ud = td>=1;
ud = td*0;
sys_d1 = ss(Ad1,Bd1,Cd,Dd,T);
sys_d2 = ss(Ad2,Bd2,Cd,Dd,T);
y_d1 = lsim(sys_d1,ud,td,x0);
y_d2 = lsim(sys_d2,ud,td,x0);

% plot system responses
figure(211);
plot(t,y_c(:,1),'Color',sstgreen,'LineStyle','-');
hold on;
plot(t,y_c(:,2),'Color',sstgreen,'LineStyle','--');
stairs(td,y_d1(:,1),'Color',sstblue,'LineStyle','-');
stairs(td,y_d1(:,2),'Color',sstblue,'LineStyle','--');
stairs(td,y_d2(:,1),'Color',sstlightblue,'LineStyle','-');
stairs(td,y_d2(:,2),'Color',sstlightblue,'LineStyle','--');
hold off;
grid on;
xlabel('t [s]');
legend('p(t)','v(t)','p_{ex}(t_k)','v_{ex}(t_k)','p_{ap}(t_k)','v_{ap}(t_k)');
title('Simple car: continuous-time, exact and approx. discretization');

% stability analysis:
[V_d,Lbd_d] = eig(Ad2),
[M,J] = jordan(Ad2),
sys_d = ss(Ad2,Bd2,Cd,Dd,T);
tf_d = tf(sys_d);

td = (0:T:10)';
ud3 = td>=1;
y_d3 = lsim(sys_d,ud3,td,x0);

% plot systems
figure(232);
stairs(td,ud3,'Color',sstgreen);
hold on;
stairs(td,y_d3(:,1),'Color',sstblue,'LineStyle','-');
stairs(td,y_d3(:,2),'Color',sstlightblue,'LineStyle','--');
hold off;
grid on;
xlabel('t [s]');
legend('u(t_k)','p(t_k)','v(t_k)');
title('Simple car: step response')

figure(233);
pzplot(sys_d,'b');
grid on;
axis equal;
% change marker size
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
    set(a(i),'color',sstblue) %change marker size
end
title('Simple car: pole zero map');

% analyze controllaboloty and observability
Mc = ctrb(Ad2,Bd2),
rank(Mc),
Mo = obsv(Ad2,Cd),
rank(Mo),

%% Exercise 1.4: no solution provided.

%% Exercise 1.5: critical points of two functions
[X1,X2] = meshgrid(-4:.2:4);

f1 = 9 - 2*X1 + 4*X2 - X1.^2 - 4*X2.^2;
f2 = 2*X1.^3 + X1.*X2.^2 + 5*X1.^2 + X2.^2;

% compute the gradient:
f1x1 = -2 - 2*X1;
f1x2 =  4 - 8*X2;
xopt1 = [-1;1/2];
fopt1 = 9 - 2*xopt1(1) + 4*xopt1(2) - xopt1(1)^2 - 4*xopt1(2)^2;

f2x1 = 6*X1.^2 + 10*X1 + X2.^2;
f2x2 = 2*X1.*X2 + 2*X2;
xopt2 = [0 , -5/3 , -1 , -1
         0 ,  0   ,  2 , -2 ];
fopt2 = 2*xopt2(1,:).^3 + xopt2(1,:).*xopt2(2,:).^2 + 5*xopt2(1,:).^2 + xopt2(2,:).^2;

Hf2(:,:,1) = [10,0;0,2];
Hf2(:,:,2) = [-10,0;0,-4/3];
Hf2(:,:,3) = [-2,4;4,0];
Hf2(:,:,4) = [-2,-4;-4,0];
lbd_Hf2 = [eig(Hf2(:,:,1)),eig(Hf2(:,:,2)),eig(Hf2(:,:,3)),eig(Hf2(:,:,4))];

figure(1511);
subplot(121);
mesh(X1,X2,f1);
hold on;
plot3(xopt1(1),xopt1(2),fopt2,'r*');
hold off;
xlabel('x_1');
ylabel('x_2');
zlabel('f_1(x)');
subplot(122);
contour(X1,X2,f1,30);
hold on;
plot(xopt1(1),xopt1(2),'r*');
quiver(X1,X2,f1x1*10,f1x2*10,'k');
hold off;
axis equal;
grid on;
legend('f_1(x)','x^*','f_{1_x}');
xlabel('x_1');
ylabel('x_2');

figure(1512);
subplot(121);
mesh(X1,X2,f2);
hold on;
plot3(xopt2(1,:),xopt2(2,:),fopt2,'r*');
hold off;
xlabel('x_1');
ylabel('x_2');
zlabel('f_1(x)');
subplot(122);
contour(X1,X2,f2,30);
hold on;
plot(xopt2(1,:),xopt2(2,:),'r*');
quiver(X1,X2,f2x1*10,f2x2*10,'k');
hold off;
axis equal;
grid on;
legend('f_2(u)','x^*','f_{2_x}');
xlabel('x_1');
ylabel('x_2');

%% Exercise 1.6: 
[U1,U2] = meshgrid(-2:.1:2);

Q = [1 1;1 2];
s = [0;1];
u_opt = -Q\s;
L_opt = -1/2*s'*(Q\s);

L = 1/2*(Q(1,1)*U1.^2 + Q(2,2)*U2.^2 + 2*Q(1,2)*U1.*U2) + s(1)*U1 + s(2)*U2;

% compute the gradient:
Lu1 = Q(1,1)*U1 + Q(1,2)*U2 + s(1);
Lu2 = Q(2,2)*U2 + Q(1,2)*U1 + s(2);

% using optimization solver:
Lfun = @(u) 1/2*(Q(1,1)*u(1)^2 + Q(2,2)*u(2)^2 + 2*Q(1,2)*u(1)*u(2)) + s(1)*u(1) + s(2)*u(2);
u0 = [1,1];
[uopt,Lopt,exitflag,output,grad,hessian] = fminunc(Lfun,u0);

figure(1611);
subplot(121);
mesh(U1,U2,L);
hold on;
plot3(u_opt(1),u_opt(2),L_opt,'r*');
hold off;
xlabel('u_1');
ylabel('u_2');
zlabel('L(u)');
subplot(122);
contour(U1,U2,L,50);
hold on;
plot(u_opt(1),u_opt(2),'r*');
quiver(U1,U2,Lu1*10,Lu2*10,'k');
hold off;
axis equal;
grid on;
legend('L(u)','u^*','L_u');
xlabel('u_1');
ylabel('u_2');

%% Exercise 1.7: no solution provided.

%% Exercise 1.8: iterative optimization algorithm

x = [0;0];
f_prev = inf;
fprintf(' \ni , x^i         , f^i    , ||Df|| , |f^i-f^(i-1)| , alpha \n');
for i = 1:11
    f = (x(1)-3)^4 + (x(1) -3*x(2))^2;
    fx = [   4*(x(1) - 3)^3 + 2*(x(1) - 3*x(2))
            -6*(x(1) - 3*x(2)) ];
    df = abs(f-f_prev);
    alpha = 0.9^i;
    f_prev = f;    
    fprintf('%i , [%.2f , %.2f] , %.3f , %.2f , %.2f , %.2f \n',i-1,x(1),x(2),f,norm(fx),df,alpha);
    
    x = x - alpha*fx/norm(fx);
end

%% Exercise 1.9: no solution provided.

%% Exercise 1.10: constrained optimization example

[X,U] = meshgrid(-2:.1:2);

Q = [1 1;1 2];
s = [0;1];
L = 1/2*(Q(1,1)*X.^2 + Q(2,2)*U.^2 + 2*Q(1,2)*X.*U) + s(1)*X + s(2)*U;
Lx = Q(1,1)*X + Q(1,2)*U + s(1);
Lu = Q(2,2)*U + Q(1,2)*X + s(2);

% plot equality constraint:
Xf = -2:0.1:2;
Uf = Xf-3;
idxs = find(Uf <= 2 & Uf>=-2);

% compute optimal solution from 1st order conditions of optimality
% x -  u           =  3
% x +  u + \lambda =  0
% x + 2u - \lambda = -1
A = [1 -1 0;1 1 1;1 2 -1]; b = [3;0;-1];
xul = A\b;
xopt = xul(1);
uopt = xul(2);

xu_opt_unc = -Q\s;

figure(2314);
contour(X,U,L);
hold on;
plot(Xf(idxs),Uf(idxs),'-.','Linewidth',2,'Color',sstlightblue);
plot(xopt,uopt,'r*');
plot(xu_opt_unc(1),xu_opt_unc(2),'s','Color',sstgreen);
quiver(X,U,Lx*10,Lu*10,'k');
hold off;
axis equal;
grid on;
xlabel('x');
ylabel('u');
legend('L(x,u)','f(x,u)=0','[x^*,u^*]','[x^*,u^*]_{unc}','[L_x,L_u]');

%% Exercise 1.11: no solution provided.

%% Exercise 1.12: optimal control example

N = 3;
P(N) = 8;
Q = 1;
R = 0.5;
A = 0.5;
B = 2;
x(1) = 1;

Kminus(N) = (B'*P(N)*B+R)^(-1)*B'*P(N)*A;
fprintf('\nk , P(k)    , K(k-1) \n');
fprintf('%i , %.5f , %.5f \n',N,P(N),Kminus(N));
for k = N-1:-1:1
    P(k) = A'*(P(k+1) - P(k+1)*B*(B'*P(k+1)*B+R)^(-1)*B'*P(k+1))*A + Q;
    Kminus(k) = (B'*P(k)*B+R)^(-1)*B'*P(k)*A;
    fprintf('%i , %.5f , %.5f \n',k,P(k),Kminus(k));
end

fprintf('\nk , x(k) , u(k-1) \n');
for k = 2:N+1
    uminus(k-1) = -Kminus(k-1)*x(k-1);
    fprintf('%i , %.5f , %.5f \n',k-2,x(k-1),uminus(k-1));
    x(k) = (A-B*Kminus(k-1))*x(k-1);
end

%% Exercise 1.13: no solution provided.


