%% Simulation

close all;
nk = 600;%250;
TU = 1:nk;
TX = 1:nk+1;
%TX = 0:T:(nk+1)*T - T;
Tref = 1:nk+N;
ref = 10 * square(0.0002*Tref, 0.79);
%ref = [ref; ref];

%player1
x01 = [xd0(1:2)*0 ; Ci{1}*xd0(1:2)];
U1 = zeros(nu1,N,nk);
X1(:,1) = x01;

%player2
x02 = [xd0(2:4)*0 ; Ci{2}*xd0(2:4)];
U2 = zeros(nu2,N,nk);
X2(:,1) = x02;
% 
% %player3
% x03 = [xd0(4:6)*0 ; Ci{3}*xd0(4:6)];
% U3 = zeros(nu2,N,nk);
% X3(:,1) = x03;

% %player4
% x04 = [xd0(6:8)*0 ; Ci{4}*xd0(6:8)];
% U4 = zeros(nu2,N,nk);
% X4(:,1) = x04;
% 
% %player5
% x05 = [xd0(8:10)*0 ; Ci{5}*xd0(8:10)];
% U5 = zeros(nu2,N,nk);
% X5(:,1) = x05;
% 
% %player6
% x06 = [xd0(10:12)*0 ; Ci{6}*xd0(10:12)];
% U6 = zeros(nu2,N,nk);
% X6(:,1) = x06;
% 
% %player7
% x07 = [xd0(12:14)*0 ; Ci{7}*xd0(12:14)];
% U7 = zeros(nu2,N,nk);
% X7(:,1) = x07;
% 
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
u(:,1) = zeros(2,1);

for k = 2:nk
    Yb = ref(:,k:k+N)';

    %Player 1
    Dxdk1 = Xd(1:2,k)-Xd(1:2,k-1);
    X1(:,k) = [ Dxdk1; Ci{1}*Xd(1:2,k)];
    x1k = X1(:,k);
    
    %Player 2
    Dxdk2 = Xd(2:4,k)-Xd(2:4,k-1);
    X2(:,k) = [ Dxdk2; Ci{2}*Xd(2:4,k)];
    x2k = X2(:,k);
    
%     %Player 3
%     Dxdk3 = Xd(4:6,k)-Xd(4:6,k-1);
%     X3(:,k) = [ Dxdk3; Ci{3}*Xd(4:6,k)];
%     x3k = X3(:,k);
    
%     %Player 4
%     Dxdk4 = [Xd(6,k);Xd(7,k);Xd(8,k)]-[Xd(6,k-1);Xd(7,k-1);Xd(8,k-1)];
%     X4(:,k) = [ Dxdk4; Ci{4}*[Xd(6,k);Xd(7,k);Xd(8,k)]];
%     x4k = X4(:,k);
%     
%     %Player 5
%     Dxdk5 = [Xd(8,k);Xd(9,k);Xd(10,k)]-[Xd(8,k-1);Xd(9,k-1);Xd(10,k-1)];
%     X5(:,k) = [ Dxdk5; Ci{5}*[Xd(8,k);Xd(9,k);Xd(10,k)]];
%     x5k = X5(:,k);    
%     
%     %Player 6
%     Dxdk6 = [Xd(10,k);Xd(11,k);Xd(12,k)]-[Xd(10,k-1);Xd(11,k-1);Xd(12,k-1)];
%     X6(:,k) = [ Dxdk6; Ci{6}*[Xd(10,k);Xd(11,k);Xd(12,k)]];
%     x6k = X6(:,k);     
%     
%     %Player 7
%     Dxdk7 = [Xd(12,k);Xd(13,k);Xd(14,k)]-[Xd(12,k-1);Xd(13,k-1);Xd(14,k-1)];
%     X7(:,k) = [ Dxdk7; Ci{7}*[Xd(12,k);Xd(13,k);Xd(14,k)]];
%     x7k = X7(:,k);    
%     
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
    

    
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
%     U3p = reshape( U3(:,:,k-1) ,[],1);
%     U4p = reshape( U4(:,:,k-1) ,[],1);
%     U5p = reshape( U5(:,:,k-1) ,[],1);
%     U6p = reshape( U6(:,:,k-1) ,[],1);
%     U7p = reshape( U7(:,:,k-1) ,[],1);
%     U8p = reshape( U8(:,:,k-1) ,[],1);
%     U9p = reshape( U9(:,:,k-1) ,[],1);
%     U10p = reshape( U10(:,:,k-1) ,[],1);
    
    %wu1 = [-U_min1 + M41*U1p;U_max1 - M41*U1p];
    %wu2 = [-U_min2 + M42*U2p;U_max2 - M42*U2p];
    %wu1 = [-U_min1 + M41*u(1,k-1);U_max1 - M41*u(1,k-1)];
    %wu2 = [-U_min2 + M42*u(2,k-1);U_max2 - M41*u(2,k-1)];
    %wy1 = [-Y_min1 + Fb1*x1k; Y_max1 - Fb1*x1k];
    %wy2 = [-Y_min2 + Fb2*x2k; Y_max2 - Fb2*x2k];
    %wr1 = [wdu1;wu1];
    %wr2 = [wdu2;wu2];

    % Get optimal sequence for player 1
    St1 = Sx{1,1}*x1k - Sy{1,1}*Yb + Sx{1,2}*x2k - Sy{1,2}*Yb + Su{1,2}*U2p;
    [U1o,J1o,exitflag] = quadprog(Rt{1},St1,[],[]);
    if exitflag<0
        error('Problems in the Optimization problem (player 1).');
    end
    
    % Get optimal sequence for player 2
    St2 = Sx{2,1}*x1k - Sy{2,1}*Yb + Sx{2,2}*x2k - Sy{2,2}*Yb + Su{2,1}*U1p ;
    [U2o,J2o,exitflag] = quadprog(Rt{2},St2,[],[]);
    if exitflag<0
        error('Problems in the Optimization problem (player 2).');
    end
    
%     % Get optimal sequence for player 3
%     St3 = Sx{3,2}*x2k - Sy{3,2}*Yb + Sx{3,3}*x3k - Sy{3,3}*Yb + Su{3,2}*U2p ;
%     [U3o,J3o,exitflag] = quadprog(Rt{3},St3,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 3).');
%     end
    
% %     % Get optimal sequence for player 4
%     St4 = Sx{4,3}*x3k - Sy{4,3}*Yb + Sx{4,4}*x4k - Sy{4,4}*Yb + Su{4,3}*U3p ;
%     [U4o,J4o,exitflag] = quadprog(Rt{4},St4,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 4).');
%     end
%     
%     % Get optimal sequence for player 5
%     St5 = Sx{5,4}*x4k - Sy{5,4}*Yb + Sx{5,5}*x5k - Sy{5,5}*Yb + Su{5,4}*U4p ;
%     [U5o,J5o,exitflag] = quadprog(Rt{5},St5,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 5).');
%     end
%     
%     % Get optimal sequence for player 6
%     St6 = Sx{6,5}*x5k - Sy{6,5}*Yb + Sx{6,6}*x6k - Sy{6,6}*Yb + Su{6,5}*U5p ;
%     [U6o,J6o,exitflag] = quadprog(Rt{6},St6,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 6).');
%     end
%     
%     % Get optimal sequence for player 7
%     St7 = Sx{7,6}*x6k - Sy{7,6}*Yb + Sx{7,7}*x7k - Sy{7,7}*Yb + Su{7,6}*U6p ;
%     [U7o,J7o,exitflag] = quadprog(Rt{7},St7,[],[]);
%     if exitflag<0
%         error('Problems in the Optimization problem (player 7).');
%     end   
%     
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
%         error('Problems in the Optimization problem (player 10).');
%     end    
    
    
    
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
%     U3pp = U3p;
%     U4pp = U4p;
%     U5pp = U5p;
%     U6pp = U6p;
%     U7pp = U7p;
%     U8pp = U8p;
%     U9pp = U9p;
%     U10pp = U10p;
    
    
    for p = 1:np
        U1pp = w(1)*U1o + (1-w(1))*U1pp;
        U2pp = w(2)*U2o + (1-w(2))*U2pp;
%         U3pp = w(3)*U3o + (1-w(3))*U3pp;
%         U4pp = w(4)*U4o + (1-w(4))*U4pp;
%         U5pp = w(5)*U5o + (1-w(5))*U5pp;
%         U6pp = w(6)*U6o + (1-w(6))*U6pp;
%         U7pp = w(7)*U7o + (1-w(7))*U7pp;
%         U8pp = w(8)*U8o + (1-w(8))*U3pp;
%         U9pp = w(8)*U9o + (1-w(9))*U9pp;
%         U10pp = w(10)*U10o + (1-w(10))*U10pp;
    end
    U1(:,:,k) = reshape( U1pp ,nu1,N);
    U2(:,:,k) = reshape( U2pp ,nu2,N);
%     U3(:,:,k) = reshape( U3pp ,nu2,N);
%     U4(:,:,k) = reshape( U4pp ,nu2,N);
%     U5(:,:,k) = reshape( U5pp ,nu2,N);
%     U6(:,:,k) = reshape( U6pp ,nu2,N);
%     U7(:,:,k) = reshape( U7pp ,nu2,N);
%     U8(:,:,k) = reshape( U8pp ,nu2,N);
%     U9(:,:,k) = reshape( U9pp ,nu2,N);
%     U10(:,:,k) = reshape( U10pp ,nu2,N);
    
    % apply first value at each player
    u1k = U1(:,1,k);
    u2k = U2(:,1,k);
%     u3k = U3(:,1,k);
%     u4k = U4(:,1,k);
%     u5k = U5(:,1,k);
%     u6k = U6(:,1,k);
%     u7k = U7(:,1,k);
%     u8k = U8(:,1,k);
%     u9k = U9(:,1,k);
%     u10k = U10(:,1,k);
    u(:,k) = [u1k+u(1,k-1)
              u2k+u(2,k-1) 
%               u3k+u(3,k-1)
%               u4k+u(4,k-1)
%               u5k+u(5,k-1)
%               u6k+u(6,k-1)
%               u7k+u(7,k-1)
%               u8k+u(8,k-1)
%               u9k+u(9,k-1)
%               u10k+u(10,k-1)
                             ];
    
    % simulate system for distributed MPC
    Xd(:,k+1) = A*Xd(:,k) + B*u(:,k); %simulate joint system
end

%%
figure(7101);


grid on;
hold on;
plot(Xd(1,:),Xd(2,:),'-','Color',sstblue);
plot(Xd(3,:), Xd(4,:), '-','Color',sstgray);
% plot(Xd(5,:), Xd(6,:), '-','Color',sstgreen);
% plot(Xd(7,:), Xd(8,:), '-','Color','red');
% plot(Xd(9,:), Xd(10,:), '-','Color','magenta');

hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('car 1','car 2', 'car 3', 'car 4', 'car 5', 'car 6', 'car 7', 'car 8', 'car 9', 'car 10');
title('Phase plot');

figure(7102);
plot(Tref,ref(1,:),'k+-');
grid on;
hold on;
plot(TX,Xd(1,:),'-','Color',sstblue);
plot(TX,Xd(2,:),'-','Color',sstgray);
plot(TX,Xd(3,:),'-','Color',sstgreen);
plot(TX,Xd(4,:),'-','Color','magenta');
% plot(TX,Xd(5,:),'-','Color','red');
% plot(TX,Xd(6,:),'-','Color','cyan');
% plot(TX,Xd(7,:),'-','Color','black');
% plot(TX,Xd(8,:),'-','Color','yellow');
% plot(TX,Xd(9,:),'-','Color','green');
% plot(TX,Xd(10,:),'-','Color','blue');


hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('$$r$$','car1 $$x_1$$ ','car1 $$x_2$$','car2$$x_1$$','car2 $$x_2$$', 'car3 $$x_1$$ ','car3 $$x_2$$','car4$$x_1$$','car4 $$x_2$$', 'car5 $$x_1$$ ','car5 $$x_2$$','car 6$$x_1$$','car6 $$x_2$$', 'car7 $$x_1$$ ','car7 $$x_2$$','car8$$x_1$$','car8 $$x_2$$', 'car9 $$x_1$$ ','car9 $$x_2$$','car10$$x_1$$','car10 $$x_2$$','Location','SouthEast');
title('State evolution');

figure(7103);
plot(TU,u(1,:),'-','Color',sstblue);
grid on;
hold on;
plot(TU,u(2,:),'-','Color',sstgray);
% plot(TU,u(3,:),'-','Color',sstgreen);
% plot(TU,u(4,:),'-','Color','magenta');
% plot(TU,u(5,:),'-','Color','red');
% plot(TU,u(6,:),'-','Color','cyan');
% plot(TU,u(7,:),'-','Color','black');
% plot(TU,u(8,:),'-','Color','yellow');
% plot(TU,u(9,:),'-','Color','green');
% plot(TU,u(10,:),'-','Color','blue');

hold off;


xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('$$u_1$$','$$u_2$$','$$u_3$$','$$u_4$$','$$u_5$$','$$u_6$$','$$u_7$$','$$u_8$$','$$u_9$$','$$u_{10}$$');
title('Input');


