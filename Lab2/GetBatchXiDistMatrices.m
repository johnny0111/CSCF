function [Fb,Gb11,Gb12,Qb,Rb] = GetBatchXiMatrices(A,B11,B12,C,N,P,Q,R,alpha)
% Returns matrices for batch computation of state sequence and cost 
% functional using incremental augmented representation

nxd = size(Bd,1);
nu = size(Bd,2);
ny = size(Cd,1);

A = [Ad , zeros(nxd,ny) ; Cd*Ad , eye(ny)];
B = [Bd ; Cd*Bd];
C = [zeros(ny,nxd),eye(ny)];
nx = size(B,1);

Fd = [];
Gd11 = [];
Gd12 = [];
Hd = [];
F = [];
G11 = [];
G12 = [];
H = [];
Qb = [];
Rb = [];
for n = 0:N
    % state matrices
    F = [F ; A^n];
    Fd = [Fd ; Ad^n];
    Gi1 = [];
    Gi2 = [];
    Gdi1 = [];
    Gdi2 = [];
    for m = 0:N-1
        ni = n-m-1;
        Gaux1 = A^ni*B11;
        Gaux2 = A^ni*B12;
        Gdaux1 = Ad^ni*Bd11;
        Gdaux2 = Ad^ni*Bd12;
        if ni < 0
            Gaux1 = Gaux1*0;
            Gaux2 = Gaux2*0;
            Gdaux1 = Gdaux1*0;
            Gdaux1 = Gdaux1*0;
        end
        Gi1 = [Gi1 , Gaux1];
        Gi2 = [Gi2 , Gaux2];
        Gdi1 = [Gdi1 , Gdaux1];
        Gdi2 = [Gdi2 , Gdaux2];
    end
    G11 = [G11 ; Gi1];
    G12 = [G12 ; Gi2];
    Gd11 = [Gd1 ; Gdi1];
    Gd12 = [Gd2 ; Gdi2];
    H = blkdiag(H,C);
    Hd = blkdiag(Hd,Cd);
    
    % cost matrices
    if n < N
        Qb = blkdiag(Qb,Q);
        Rb = blkdiag(Rb,R);
    else
        Qb = blkdiag(Qb,P);
    end
    
    % compute output matrices
    Fb = H*F;
    Gb11 = H*G11;
    Gb12 = H*G12;
    
    % só para não confundir
    Qb = alpha*Qb;
    Rb = alpha*Rb;    

end



