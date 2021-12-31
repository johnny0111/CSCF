function [Kn,isfinite,Nb] = GetKn(O,sys,nsteps)
% compute the n-step or maximal controllable set for system sys, given a 
% target set O in X, and the input constraints U.

    %initialization
    isfinite = 0;
    Nb = 0;
    X = sys.x.boundsToPolyhedron();
    U = sys.u.boundsToPolyhedron();
    K_prev = O;
    K = O;
    for i = 1:nsteps
        PreK = sys.reachableSet('X',K,'U',U,'direction','backward');
        K = PreK.intersect(X).minHRep();
        if K == K_prev
            isfinite = 1;
            Nb = i;
            break;
        end
        K_prev = K;
    end
    Kn = K;
end


