% Lwh 20220401
% 输入：normal_edge （数据上边界）
% 输出：
% ytick_this： 适合数据的ytick
% ylimit_max: 适合的上限，暂时仅考虑[0 max_data]

function [ylimit_max, ytick_this] = properYtick(data)

 % input top bound or input raw data
 normal_edge = max(abs(data));
 
if normal_edge>25 % 较大的y轴,十位取3,4,5的倍数
    ylimit_max = near_multiple(normal_edge/10,[3 4 5])*10;
    if mod(ylimit_max/10,4)==0 % 优先分4个tick
        ytick_this = [0:ylimit_max/4:ylimit_max];
    elseif mod(ylimit_max/10,3)==0 % 其次分3个tick
        ytick_this = [0:ylimit_max/3:ylimit_max];
    elseif mod(ylimit_max/10,5)==0 % 其次分5个tick
        ytick_this = [0:ylimit_max/5:ylimit_max];
    end
elseif normal_edge>3
    % set the yticklabel: 4/5倍数
    ylimit_max = near_multiple(normal_edge,[4 5]);
    if mod(ylimit_max,4)==0 % 优先分4个tick
        ytick_this = [0:ylimit_max/4:ylimit_max];
    elseif mod(ylimit_max,5 )==0 % 其次分5个tick
        ytick_this = [0:ylimit_max/5:ylimit_max];
    end
else
    ylimit_max = ceil(normal_edge*10)/10; % 取0.9,0.8,0.7...等
    ytick_this = [0:0.1:ylimit_max];
end


end
