% Phase trellis for the binary CPFM-FRR scheme using weakly orthogonal signal
function [P_State,P_Ip,Ga_Inx]= Get_Trellis_manual()
P_State = [4,2;1,3;2,4;3,1]; % previous state
P_Ip = [1,2;1,2;1,2;1,2]; % previous ip
Ga_Inx = [4,8;1,7;2,6;3,5]; % branch indices
end