function [  Rt1,S11x,S11y,S12x,S12y,S12u,...
            Rt2,S22x,S22y,S21x,S21y,S21u] = GetDistMPCbatch(Fb1,Gb11,Gb12,Qb1,Rb1,...
                                                            Fb2,Gb21,Gb22,Qb2,Rb2)
% compute MPC cost matrices for two player game

% MPC 1 gains
Rt1 = Rb1 + Gb11'*Qb1*Gb11 + Gb21'*Qb2*Gb21;
S11x = Gb11'*Qb1*Fb1;
S11y = Gb11'*Qb1;
S12x = Gb21'*Qb2*Fb2;
S12y = Gb21'*Qb2;
S12u = Gb11'*Qb1*Gb12 + Gb21'*Qb2*Gb22;
K11 = -inv(Rt1)*S11x;
K12 = -inv(Rt1)*S12x;
K11y = inv(Rt1)*S11y;
K12y = inv(Rt1)*S12y;
L1 = -inv(Rt1)*S12u;

% MPC 2 gains
Rt2 = Rb2 + Gb22'*Qb2*Gb22 + Gb12'*Qb1*Gb12;
S22x = Gb22'*Qb2*Fb2;
S22y = Gb22'*Qb2;
S21x = Gb12'*Qb1*Fb1;
S21y = Gb12'*Qb1;
S21u = Gb22'*Qb2*Gb21 + Gb12'*Qb1*Gb11;
K22 = -inv(Rt2)*S22x;
K21 = -inv(Rt2)*S21x;
K22y = inv(Rt2)*S22y;
K21y = inv(Rt2)*S21y;
L2 = -inv(Rt2)*S21u;

end

