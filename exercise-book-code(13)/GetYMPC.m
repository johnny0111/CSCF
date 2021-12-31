function [Rt,Sx,Sy,Su,K,Ky,L] = GetYMPC(Fb,Gb,Qb,Rb)
% compute MPC cost matrices for m-player game (using cell arrays/matrices)

    m = length(Fb);
    for i = 1:m
        Rt{i} = Rb{i};
        for j = 1:m
            Rt{i} = Rt{i} + Gb{j,i}'*Qb{j}*Gb{j,i};
            Sy{i,j} = Gb{j,i}'*Qb{j};
            Sx{i,j} = Gb{j,i}'*Qb{j}*Fb{j};
            Su{i,j} = Gb{i,i}'*Qb{i}*Gb{i,j} + Gb{j,i}'*Qb{j}*Gb{j,j};
        end
        % compute gains also
        for j = 1:m
            K{i,j} = -Rt{i}^(-1)*Sx{i,j};
            Ky{i,j} = Rt{i}^(-1)*Sy{i,j};
            L{i,j} = -Rt{i}^(-1)*Su{i,j};
        end
    end

    if m==1
        Rt = Rt{1};
        Sx = Sx{1};
        Sy = Sy{1};
        Su = Su{1};
        K  = K{1};
        Ky = Ky{1};
        L  = L{1};
    end
end

