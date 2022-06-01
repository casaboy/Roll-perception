% Lwh 20220401
% use to find the bin size and bin edge
% adjust the normal binedge to "整数”， exclude outlier
% 为了画图好看，规则：
% 1. >24: 取30，40，50等
% 2. >2: 取最接近的3/4的倍数, 例如2.1=3, 2.5=3, 5=6, 7=8
% 3. >=1 & <=2: 取[1, 1.5, 2]，例如1=1, 1.1=1.5, 1.6=2
% 4. <1: 取0.9,0.8,0.7...等
%
% input
% data: data used for plot hist
% setsize: desired each bin size, only for "small" or "large"
%
% output
% binedge： 适合数据的左右边界 （normal_edge为右边界）
% bin_size: 每个bin的大小
% binmid: 每个bin的中心值
% xtick_this： 适合数据的xtick



% ***  此文件参数只根据我自己的数据设定的规则，可能不具有广泛适用能力 ***

function [binedge, binmid, xtick_this, bin_size] = properHistBin(data,setsize)

if nargin == 1
    setsize = 'small';
end

value_max = max(abs(data));

%% 1. find normal_edge (左右边界，取正负对称值)
if value_max > 24 % 较大的x轴
    normal_edge = ceil(value_max/10)*10; % 取30，40，50等
    
elseif  value_max > 2 % 取最接近的3/4的倍数；目标：使得xtick包括0一共3或4个数字，且每个bin的宽度比较易读
    normal_edge = near_multiple(value_max,[3 4]);
    
elseif  value_max >=1 && value_max <=2 % use 1, 1.5
    set_max = [1 1.5 2];
    temp = set_max - value_max;
    normal_edge = set_max(temp == min(temp(temp>=0))); % 最接近set_max且小于等于set_max的值
    
else %<1, for motion type task, typical PSE is very small
    normal_edge = ceil(value_max*10)/10; % 取0.9,0.8,0.7...等
    
end

%% 2. 使每个bin的宽度比较易读
if strcmpi(setsize,'small')
    if normal_edge==6 || normal_edge==9 || normal_edge==15 % 12 暂时用10个bin
        bin_size = normal_edge/9;
    else
        bin_size = normal_edge/10;
    end
    
elseif strcmpi(setsize,'large')
%      if normal_edge==6 || normal_edge==9 || normal_edge==15 % 12 暂时用10个bin
%         bin_size = normal_edge/3;
%     else
        bin_size = normal_edge/5;
%      end
else
    disp('就你逼事多');
    keyboard 
end

%% 3. binedge
binedge = -normal_edge: bin_size: normal_edge; % normal data edge
binmid = binedge(1:end-1) + bin_size/2; % for bar plot

%% 4. xticklabel
xtick_this = properXtick(normal_edge);

end