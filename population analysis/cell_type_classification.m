% cell type classification 2020 Lwh
% use in Plot_Spiral_Tuning_Lwh2.m and save the cell type

function [cell_type,mt_type,num] = cell_type_classification(data,select_all,p_cri)
num1 = 0; % only heading response significant higher than spon and with heading tuning (Heading response and heading tuning)
num2 = 0; % only heading response significant higher than spon and without heading tuning (Heading response no heading tuning)
num3 = 0; % only rotation response significant higher than spon and with rotation tuning (Rotation response and rotation tuning)
num4 = 0; % only rotation response significant higher than spon and without rotation tuning (Rotation response no rotation tuning)
num5 = 0; % both response significant higher than spon and only with heading tuning (Both response and heading tuning)
num6 = 0; % both response significant higher than spon and only with rotation tuning (Both response and rotation tuning)
num7 = 0; % both response significant higher than spon and with both tuning (Both response and both tuning)
num8 = 0; % both response significant higher than spon but without any tuning (Both response no tuning)
num9 = 0; % both response not significant higher than spon (No response)
num10 = 0; % one response significant higher than spon and with motion type tuning
num11 = 0; % one response significant higher than spon and without motion type tuning

N = length(data);
cell_type = [];

if isempty(select_all)
    select_all = true(length(data),1);
end

% 1. 先判断heading/rotation response是否高于spon
% 2. 再判断heading tuning/rotation tuning是否显著
for n = 1:N
    if select_all(n) == true
        if data(n).p_stim < 0.05  % h or r显著高于spon
            % mt_type: heading plane 和 rotation plane有显著差异
            if data(n).p_turn_mt < 0.05
                mt_type(n) = 1;
                num10 = num10 + 1; % one response significant higher than spon and with motion type tuning
            else
                mt_type(n) = 0;
                num11 = num11 + 1; % one response significant higher than spon and without motion type tuning
            end
        elseif data(n).p_stim >= 0.05 % not above spon
            mt_type(n) = 0;
        else
            mt_type(n) = nan;
        end
        
        %  分别判断h or r是否显著高于spon
        if data(n).p_stim_h < 0.05 && data(n).p_stim_r >= 0.05 % 只有h显著高于spon，此时不考虑rotation
            if data(n).p_turn_h < p_cri % h有tuning
                cell_type(n) = 1; % translation cell
                num1 = num1 + 1; % Heading response and heading tuning
            else
                cell_type(n) = 4; % no tuning but above spon
                num2 = num2 + 1; % Heading response no heading tuning
            end
        elseif data(n).p_stim_h >= 0.05 && data(n).p_stim_r < 0.05 % 只有r显著高于spon，此时不考虑heading
            if data(n).p_turn_r < p_cri % r有tuning
                cell_type(n) = 2; % rotation cell
                num3 = num3 + 1; % Rotation response and rotation tuning
            else
                cell_type(n) = 4; % no tuning but above spon
                num4 = num4 + 1; % Rotation response no rotation tuning
            end
        elseif data(n).p_stim_h < 0.05 && data(n).p_stim_r < 0.05 % h，r同时显著高于spon
            if data(n).p_turn_h < p_cri && data(n).p_turn_r >= p_cri % 只有heading tuning
                cell_type(n) = 1; % translation cell
                num5 = num5 + 1; % Both response and heading tuning
            elseif data(n).p_turn_h >= p_cri && data(n).p_turn_r < p_cri % 只有rotation tuning
                cell_type(n) = 2; % rotation cell
                num6 = num6 + 1; % Both response and rotation tuning
            elseif data(n).p_turn_h < p_cri && data(n).p_turn_r < p_cri % 同时两个tuning
                cell_type(n) = 3; % spiral cell
                num7 = num7 + 1; % Both response and both tuning
            else
                cell_type(n) = 4; % no tuning but above spon
                num8 = num8 + 1; % Both response no tuning
            end
        elseif data(n).p_stim_h >= 0.05 && data(n).p_stim_r >= 0.05 %
            cell_type(n) = 5; % not above spon (non-excitory)
            num9 = num9 + 1; % No response
        else
            cell_type(n) = nan; % 没测tuning
        end
    end
end

cell_type = cell_type';
mt_type = mt_type';
num = [num1,num2,num3,num4,num5,num6,num7,num8,num9,num10,num11];
end