function [X0,isfinite,Nb,Xi] = GetXi(Xf,sys,N)
% compute the N-step feasibility set for system sys, given the final 
% invariant constraint set Xf, state constraint X, and 
% input constraint U (the last two given inside sys).

    isfinite = 0;
    Nb = 0;
    X = sys.x.boundsToPolyhedron();
    U = sys.u.boundsToPolyhedron();
    Xi = Polyhedron;
    Xi(N+1) = Xf;
    i = N;
    while i >= 1
        PreX = sys.reachableSet('X',Xi(i+1),'U',U,'direction','backward');
        Xi(i) = PreX.intersect(X).minHRep();
        if Xi(i) == Xi(i+1)
            isfinite = 1;
            Nb = i;
            break;
        end
        i = i-1;
    end
    X0 = Xi(i+1);
end