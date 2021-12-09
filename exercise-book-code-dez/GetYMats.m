function [Fb,Gb,Qb,Rb,F,G,H] = GetYMats(A,B,C,N,P,Q,R,alpha)
% Returns matrices for batch computation of state sequence and cost 
% functional for distributed m-player problems (using cell arrays/matrices)

    if ~exist('alpha','var')
        alpha = 1;
    end

    if ~iscell(A)
        A = {A};
        B = {B};
        C = {C};
    end
    m = length(A);

    for i = 1:m

        F{i} = [];
        H{i} = [];
        Qb{i} = [];
        Rb{i} = [];
        for n = 0:N
            % state matrices
            F{i} = [F{i} ; A{i}^n];
            H{i} = blkdiag(H{i},C{i});

            % cost matrices
            if exist('P','var') && exist('Q','var') && exist('R','var')
                if n < N
                    Qb{i} = blkdiag(Qb{i},Q);
                    Rb{i} = blkdiag(Rb{i},R);
                else
                    Qb{i} = blkdiag(Qb{i},P);
                end
            end
        end
        Qb{i} = alpha*Qb{i};
        Rb{i} = alpha*Rb{i};
        Fb{i} = H{i}*F{i};
        
        for j = 1:m
            G{i,j} = [];
            for n = 0:N
                Gnj = [];
                for mu = 0:N-1
                    nm = n-mu-1;
                    Gauxj = A{i}^nm*B{i,j};
                    if nm < 0
                        Gauxj = Gauxj*0;
                    end
                    Gnj = [Gnj , Gauxj];
                end
                G{i,j} = [G{i,j} ; Gnj];
            end
            Gb{i,j} = H{i}*G{i,j};
        end
        
    end

end

