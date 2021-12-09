function [ Fb,Gb1,Gb2,Qb,Rb] = GetBatchYdistMatrices(A,B1,B2,C,N,P,Q,R,alpha)
                                                        
    F = [];
    G1 = [];
    G2 = [];
    H = [];
    Qb = [];
    Rb = [];
    for n = 0:N
        % state matrices
        F = [F ; A^n];
        Gi1 = [];
        Gi2 = [];
        for m = 0:N-1
            ni = n-m-1;
            Gaux1 = A^ni*B1;
            Gaux2 = A^ni*B2;
            if ni < 0
                Gaux1 = Gaux1*0;
                Gaux2 = Gaux2*0;
            end
            Gi1 = [Gi1 , Gaux1];
            Gi2 = [Gi2 , Gaux2];
        end
        G1 = [G1 ; Gi1];
        G2 = [G2 ; Gi2];
        H = blkdiag(H,C);

        % cost matrices
        if exist('P','var') && exist('Q','var') && exist('R','var')
            if n < N
                Qb = blkdiag(Qb,Q);
                Rb = blkdiag(Rb,R);
            else
                Qb = blkdiag(Qb,P);
            end
        end
    end

    % compute output matrices
    Fb = H*F;
    Gb1 = H*G1;
    Gb2 = H*G2;
    Qb = alpha*Qb;
    Rb = alpha*Rb;

end


