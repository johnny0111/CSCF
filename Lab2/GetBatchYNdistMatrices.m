function [Fb,Gb,Qb,Rb] = GetBatchYNdistMatrices(Ai,Bi,Ci,N,Pi,Qi,Ri,alphai)
    
    np = size(Ai,2);
    F = cell(1,np);
    G = cell(np,np);
    H = cell(1,np);
    Qb = cell(1,np);
    Rb = cell(1,np);

    
    for i = 1:np
        for n = 0:N
            % state matrices
            F{i} = [F{i} ; Ai{i}^n];
            
            for j = 1:np
                Giii = [];
                for m = 0:N-1
                    
                    ni = n-m-1;
                    Gauxii = Ai{i}^ni*Bi{i,j};
                    if ni < 0
                        Gauxii = Gauxii*0;
                    end
                    Giii = [Giii , Gauxii];
                    
                end
                G{i,j} = [G{i,j} ; Giii];
            end
            
            
            H{i} = blkdiag(H{i},Ci{i});
            
            % cost matrices
            if exist('Pi','var') && exist('Qi','var') && exist('Ri','var')
                if n < N
                    Qb{i} = alphai*blkdiag(Qb{i},Qi(i));
                    Rb{i} = alphai*blkdiag(Rb{i},Ri(i));
                else
                    Qb{i} = alphai*blkdiag(Qb{i},Pi);
                end
            end
        end
        % compute output matrices
        
        
        
    end
    for i=1:np
        Fb{i} = H{i}*F{i};
        for j = 1:np
            Gb{i,j} = H{i}*G{i,j};
        end
    end
    
    
end