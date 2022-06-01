% Lwh 20220401
% 输入：normal_edge （数据右边界）
% 输出：xtick_this： 适合数据的xtick
% 适合左右对称的轴

function xtick_this = properXtick(data)

%% find the proper right bopund 
if size(data,1) == 1 && size(data,2) == 1 % input right bound
    normal_edge = data;
    
else % input raw data
    % find the proper right bopund
    value_max = max(abs(data));
    if value_max > 24 % 较大的x轴
        normal_edge = ceil(value_max/10)*10; % 取30，40，50等
        
        % 解决50，70无法被3或4整除
        if normal_edge == 50
            normal_edge = 60;
        elseif normal_edge == 70
            normal_edge = 80;
        end
        
    elseif  value_max > 2 % 取最接近的3/4的倍数；目标：使得xtick包括0一共3或4个数字，且每个bin的宽度比较易读
        normal_edge = near_multiple(value_max,[3 4]);
        
    elseif  value_max >=1 && value_max <=2 % use 1, 1.5
        set_max = [1 1.5 2];
        temp = set_max - value_max;
        normal_edge = set_max(temp == min(temp(temp>=0))); % 最接近set_max且小于等于set_max的值
        
    else %<1, for motion type task, typical PSE is very small
        normal_edge = ceil(value_max*10)/10; % 取0.9,0.8,0.7...等
    end
end

%% find the proper xitck
if mod(normal_edge,3)==0 % 优先分3个tick
    xtick_this = [-normal_edge:normal_edge/3:normal_edge];
    
elseif mod(normal_edge,4)==0 % 其次分4个tick
    xtick_this = [-normal_edge:normal_edge/4:normal_edge];
    
elseif normal_edge == 2 || normal_edge == 1.5
    xtick_this = [-normal_edge:0.5:normal_edge];
else % <1
    xtick_this = [-normal_edge:0.1:normal_edge];
    
end

if normal_edge == 12 % 分4个tick
    xtick_this = [-normal_edge:normal_edge/4:normal_edge];
end

end