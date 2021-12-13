function [Rt,Sx,Sy,Su] = GetYMPC(Fb,Gb,Qb,Rb);
    np = size(Gb,1);
    Rt = cell(1,np);
    Sx = cell(np,np);
    Sy = cell(np,np);
    Su = cell(np,np);
    
    
    for i = 1:np
        Rt{i} = Rb{i};
        for j = 1:np
            Rt{i} = Rt{i} + Gb{j,i}'*Qb{j}*Gb{j,i};
            Sx{i,j} = Gb{j,i}'*Qb{j}*Fb{j};
            Sy{i,j} = Gb{j,i}'*Qb{j};
            if (i ~= j)
                Su{i,j} = Gb{i,i}'*Qb{i}*Gb{i,j} + Gb{j,i}'*Qb{j}*Gb{j,j};
                for k = 1:np
                    if(k ~= i & k~= j)
                    Su{i,j} = Su{i,j} + 
                    end
                end
            end
            
            
        end
    end

end

