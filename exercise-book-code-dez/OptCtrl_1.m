function [x,u,lambda,k] = OptCtrl_1(a,b,c,r,N,x0,rN)
%simulate class example of NL opt control
% auxiliary OC parameters
    gamma = b^2/r;
    Lamb = gamma*a^(2*(N-1))*(1-a^(-2*N))/(1-a^(-2));
    C = c*(1-a^N)/(1-a);

    % compute state, costate and control evolution
    x = zeros(1,N+2);
    u = zeros(1,N+2);
    lambda = zeros(1,N+2);
    x(1) = x0;
    lambda(1) = a^(N-1)/(1+Lamb)*(a^N*x0 - rN + C);
    for k = 1:N
        lambda(k+1) = a^(N-k-1)/(1+Lamb)*(a^N*x0 - rN + C);
        u(k) = -b/r*lambda(k+1);
        x(k+1) = a*x(k) - gamma*lambda(k+1) + c;
    end
    u(k+1) = u(k);
    u(k+2) = u(k);
    x(k+2) = x(k+1);
    lambda(k+2) = lambda(k+1);
    k = 0:N+1;

end

