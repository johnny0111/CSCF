function [Fbi,Gbii,Gbij,Qbi,Rbi] = GetBatchYdistMatrices(Ai,Bii,Bij,Ci,N,Pi,Qi,Ri,alphai)

    
    Fi = [];
    Gii = [];
    Gij = [];
    Hi = [];
    Qbi = [];
    Rbi = [];
    for n = 0:N
        % state matrices
        Fi = [Fi ; Ai^n];
        Giii = [];
        Giij = [];
        for m = 0:N-1
            ni = n-m-1;
            Gauxii = Ai^ni*Bii;
            Gauxij = Ai^ni*Bij;
            if ni < 0
                Gauxii = Gauxii*0;
                Gauxij = Gauxij*0;
            end
            Giii = [Giii , Gauxii];
            Giij = [Giij , Gauxij];
        end
        Gii = [Gii ; Giii];
        Gij = [Gij ; Giij];
        Hi = blkdiag(Hi,Ci);

        % cost matrices
        if exist('Pi','var') && exist('Qi','var') && exist('Ri','var')
            if n < N
                Qbi = alphai*blkdiag(Qbi,Qi);
                Rbi = alphai*blkdiag(Rbi,Ri);
            else
                Qbi = alphai*blkdiag(Qbi,Pi);
            end
        end
    end

    % compute output matrices
    Fbi = Hi*Fi;
    Gbii = Hi*Gii;
    Gbij = Hi*Gij;
    
end