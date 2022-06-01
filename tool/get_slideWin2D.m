% get slide windows data
% Lwh
% xdata,ydata is a vector
% 1D: step bin across x axis, mean y axis, use "get_slideWin"
% 2D: step bin across x and y axis, mean a square area, use for heatmap or imagesc or etc. use "get_slideWin2D"
% 2D: win = [xwin,ywin], step = [xstep ystep], xdata: 
% 边缘的数据可能被忽略！！！ 使用合适的win和step
% example for using pcolor:
% win = [2 1];
% step = [1 0.5];
% [xmid_slide, ymid_slide, z, z_sem] = get_slideWin2D(XX,YY,CC,win,step,[-10 10],[-5 5]);
% [aa,bb] = size(z);
% z_plot = z;
% z_plot(aa+1,:) = z_plot(aa,:);  % 伪造最后一行
% z_plot(:,bb+1) = z_plot(:,bb); % 伪造最后一列
% 
% x_plot = unique(xmid_slide) - step(1)/2; % 调整坐标为每个格子的左下角网格交点
% x_plot(end+1) = x_plot(end) + step(1);
% 
% y_plot = unique(ymid_slide) - step(2)/2; % 调整坐标为每个格子的左下角网格交点
% y_plot(end+1) = y_plot(end) + step(2);
% [x_plot, y_plot] = meshgrid(x_plot,y_plot);
% 
% s = pcolor(x_plot, y_plot,z_plot); % 每个格子中心对应原来xmid，ymid
% axis square
% colorbar
% s.FaceColor = 'interp'; % 跨格子插值
% SetFigure(15)

function [xmid, ymid, z, z_sem] = get_slideWin2D(xdata,ydata,zdata,win,step,xedges,yedges)

% [n,m] = size(xdata);
% turn to column
xdata = xdata(:);
ydata = ydata(:);
zdata = zdata(:);

if nargin < 6
    xedges = [min(xdata) max(xdata)]; % auto use min and max
    yedges = [min(ydata) max(ydata)]; % auto use min and max
end

if length(win)==1
    win = [win win]; % use same window in x and y
end

if length(step)==1
    step = [step step]; % use same step in x and y
end

xwind(1,1) = xedges(1);
xwind(1,2) = xedges(1)+win(1);

ywind(1,1) = yedges(1);
ywind(1,2) = yedges(1)+win(2);

if xwind(2)>xedges(2) || ywind(2)>yedges(2) % too large slide window
    keyboard
end

% get windows
i = 1;
while xwind(i, 2)<xedges(2)
    i = i + 1;
    xwind(i, 1) = xwind(i-1, 1)+step(1);
    xwind(i, 2) = xwind(i-1, 2)+step(1);
end

j = 1;
while ywind(j,2)<yedges(2)
    j = j+1;
    ywind(j,1) = ywind(j-1,1)+step(2);
    ywind(j,2) = ywind(j-1,2)+step(2);
end

if xwind(end, 2)~=xedges(2) || ywind(end, 2)~=yedges(2)
    keyboard
    disp('maybe loss edge data!');
end

for i = 1:size(xwind,1)
   for j = 1:size(ywind,1) 
       if xwind(i,2)==xedges(2)
            select = logical(xdata >= xwind(i,1) & xdata <= xwind(i,2) & ydata >= ywind(j,1) & ydata < ywind(j,2) ); % deal with x right edge
        elseif ywind(j,2)==yedges(2) 
            select = logical(xdata >= xwind(i,1) & xdata < xwind(i,2) & ydata >= ywind(j,1) & ydata <= ywind(j,2) ); % deal with y right edge
        elseif xwind(i,2)==xedges(2) && ywind(j,2)==yedges(2)
            select = logical(xdata >= xwind(i,1) & xdata <= xwind(i,2) & ydata >= ywind(j,1) & ydata <= ywind(j,2) ); % deal with x&y right edge
        else
            select = logical(xdata >= xwind(i,1) & xdata < xwind(i,2) & ydata >= ywind(j,1) & ydata < ywind(j,2) ); % normal: >= left edge, < right edge
        end
        
        xmid(j,i) = xwind(i,1) + win(1)/2;
        ymid(j,i) = ywind(j,1) + win(2)/2;
        
        if sum(select)>0
            z(j,i) = nanmean(zdata(select));
            z_sem(j,i) = nanstd(zdata(select)) / sqrt(sum(select));
        else
            z(j,i) = nan;
            z_sem(j,i) = nan;
        end
       
   end
end

end