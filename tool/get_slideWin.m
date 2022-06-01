% get slide windows data
% Lwh
% xdata,ydata is a vector
% 1D: step bin across x axis, mean y axis
% 2D: step bin across x and y axis, mean a square area, use for heatmap or
% imagesc or etc. use "get_slideWin2D"
% 边缘的数据可能被忽略！！！使用合适的win和step

function [xmid, y, y_sem, y_step_data] = get_slideWin(xdata,ydata,win,step,xedges)

if nargin == 4
    xedges = [min(xdata) max(xdata)]; % auto use min and max
end

xwind(1,1) = xedges(1);
xwind(1,2) = xedges(1)+win;

if xwind(2)>xedges(2) % too large slide window
    keyboard
end

% get windows
i = 1;
while xwind(i, 2)<xedges(2)
    i = i + 1;
    xwind(i, 1) = xwind(i-1, 1)+step;
    xwind(i, 2) = xwind(i-1, 2)+step;
end

if abs(xwind(end, 2)-xedges(2))>10e-6
    keyboard
    disp('maybe loss edge data!');
end

for i = 1:size(xwind,1)
    if xwind(i, 2)==xedges(2)
        select = logical(xdata >= xwind(i, 1) & xdata <= xwind(i, 2) ); % x <= a <= x+1
    else
        select = logical(xdata >= xwind(i, 1) & xdata < xwind(i, 2) ); % x <= a < x+1
    end
    
    xmid(i) = xwind(i, 1) + win/2; % middle of windows
    
    if sum(select)>0
        y(i) = nanmean(ydata(select));
        y_sem(i) = nanstd(ydata(select)) / sqrt(sum(select));
        
        % save y
        y_step_data_temp = ydata(select);
        y_step_data{i} = y_step_data_temp(~isnan(y_step_data_temp));
    else
        y(i) = nan;
        y_stm(i) = nan;
        y_step_data{i} = nan;
    end
end
