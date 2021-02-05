% Phase trellis for the binary CPFM-FRR scheme using weakly orthogonal signal
function [Prev_State,Prev_Ip,Outputs_prev]= Get_Trellis_manual()
Prev_State = [4,2;1,3;2,4;3,1]; % previous state
Prev_Ip = [1,2;1,2;1,2;1,2]; % previous ip
Outputs_prev = [4,8;1,7;2,6;3,5]; % branch indices
end