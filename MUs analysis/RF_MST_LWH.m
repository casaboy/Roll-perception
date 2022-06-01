% MST LWH
% 20181127
% 从Batch file中获取每个细胞的RF信息（如果有的话）
[session, hemi, xloc, yloc, depth, FILE, rf, fixation_x, fixation_y,ab_stim] = GetBatchfileInfo(1);
%[yloc, xloc, FILE, rf, eye_x, eye_y] = GetRFInfo; %得到每一个file对应的RF（有些cell可能没有测RF）   RF&File.m

% char转为num
session = cellfun(@str2num,session);
hemi = cellfun(@str2num,hemi);
xloc = cellfun(@str2num,xloc);
yloc = cellfun(@str2num,yloc);
depth = cellfun(@str2num,depth);
f_x = cellfun(@str2num,fixation_x);
f_y = cellfun(@str2num,fixation_y);

% still cell
rf = cellfun(@(x) str2num(x),rf,'UniformOutput', false);

all_num = size(session,2);

Viewing_distance=45;

if ishandle(101); close(101); end; figure(101);
set(101,'Position', [500,250, 1000,500], 'Name', 'MST RF', 'color','w');

subplot(1,2,1)
axis equal; box on; axis([-60 60 -60 60]);
line([-60 60],[0 0],'color','k','LineStyle',':'); hold on;
line([0 0],[-60 60],'color','k','LineStyle',':');
set(gca,{'xtick','ytick'},{-60:20:60,-60:20:60});
xlabel('degree');
ylabel('degree');
title('MST RF');

MST_RF_x = [];
MST_RF_y = [];
MST_RF_w = [];
MST_RF_h = [];

for i=1:all_num
    try
        MST_RF_x(i) = rf{i,1} ;
        MST_RF_y(i) = rf{i,2} ;
        MST_RF_w(i) = rf{i,3} ;
        MST_RF_h(i) = rf{i,4};
    catch
        keyboard
    end
end

true_num = 0;
for i=1:all_num
    %if MST_RF_w>50 & MST_RF_h >50
    if ~isnan(MST_RF_x(i))
    rectangle('position',[MST_RF_x(i)-MST_RF_w(i)/2 MST_RF_y(i)-MST_RF_h(i)/2, MST_RF_w(i) MST_RF_h(i)],...
        'Curvature',[0.3 0.3],'EdgeColor','k','LineWidth',1.5);
    hold on
    true_num = true_num+1;
    %end
    end
end
text(-55,55,sprintf('N=%g',true_num));


% 转化为圆的半径比较
d = [];
for i=1:all_num
    area(i) = rf{i,3} * rf{i,4};
    if ~isnan(area(i))
        d(i) = 2 * sqrt(area(i) / pi);
    else
        d(i) = NaN;
    end
end

%max and min
maxrf = find(area == nanmax(area));
minrf = find(area == nanmin(area));

rectangle('position',[MST_RF_x(maxrf)-MST_RF_w(maxrf)/2 MST_RF_y(maxrf)-MST_RF_h(maxrf)/2, MST_RF_w(maxrf) MST_RF_h(maxrf)],...
    'Curvature',[0.3 0.3],'EdgeColor','r','LineWidth',5);
hold on
rectangle('position',[MST_RF_x(minrf)-MST_RF_w(minrf)/2 MST_RF_y(minrf)-MST_RF_h(minrf)/2, MST_RF_w(minrf) MST_RF_h(minrf)],...
    'Curvature',[0.3 0.3],'EdgeColor','b','LineWidth',5);
hold off

text(-55,75,sprintf('Max:%g x %g ',MST_RF_w(maxrf),MST_RF_h(maxrf)));
text(-55,70,sprintf('Min:%g x %g ',MST_RF_w(minrf),MST_RF_h(minrf)));


% 转化为圆的半径比较
subplot(1,2,2)
xbins = [0:10:110];
h = hist(d,xbins);
hold on;
h_h = bar(xbins,h,1,'facecolor','w');

xlabel('Equivalent diameter of receptive field (degree)');
ylabel('Number of cases');
set(gca,'xlim',[0 120]);
%set(gca,'ylim',[0 20]);
%set(gca,'ytick',[0:5:20]);

% mean
mean_d = nanmean(d);
top = max(h);
plot([mean_d mean_d],[0 top+1],'k:');
plot(mean_d,top+1.5,'kv','markerfacecolor','k');
text(mean_d+4,top+1.5,num2str(mean_d),'fontsize',15);

% 包括fovea的细胞数目
% if ishandle(102); close(102); end; figure(102);
% set(102,'Position', [500,250, 700,600], 'Name', 'MST RF', 'color','w');
includfovea=0;
for i=1:all_num
    if (MST_RF_x(i)-MST_RF_w(i)/2) < -0 & (MST_RF_x(i)+MST_RF_w(i)/2) > 0 & (MST_RF_y(i)-MST_RF_h(i)/2) < -0 & (MST_RF_y(i)+MST_RF_h(i)/2) > 0
%         rectangle('position',[MST_RF_x(i)-MST_RF_w(i)/2 MST_RF_y(i)-MST_RF_h(i)/2, MST_RF_w(i) MST_RF_h(i)],...
%             'Curvature',[0.3 0.3],'EdgeColor','k','LineWidth',1.5);
%         hold on
        includfovea=includfovea+1;
    end
end

subplot(1,2,1)
text(-55,50,sprintf('Included Fovea N=%g',includfovea));
% axis equal; box on; axis([-60 60 -60 60]);
% line([-60 60],[0 0],'color','k','LineStyle',':'); hold on;
% line([0 0],[-60 60],'color','k','LineStyle',':');
% set(gca,{'xtick','ytick'},{-60:20:60,-60:20:60});
% xlabel('degree');
% ylabel('degree');
% title('Included fovea');

%% 在grid x-y平面上，看不同x时候rf的规律（从中间（sts tip处，小数值）到外（大数值））
if ishandle(102); close(102); end; figure(102);
set(102,'Position', [500,250, 1000,500], 'Name', 'MST RF', 'color','w');

% xbin = linspace(min(xloc),max(xloc),4);
% ybin = linspace(min(yloc),max(yloc),4);


xlocbin = min(xloc):2:max(xloc)+2;
ylocbin = min(yloc):3:max(yloc)+3;

% 画distribution即可
site_d = []; mean_site_d = [];
for i = 1:length(xlocbin)-1
    for j = 1:length(ylocbin)-1
        %select = find(xloc>=xlocbin(i) & xloc<xlocbin(i+1) & yloc>=ylocbin(j) & yloc<ylocbin(j+1));
        select = find((xloc==18 | xloc==19|xloc==20) );
        
        %if ~isempty(select)
        for n = 1:length(select)
            site_x = rf{select(n),1};
            site_y = rf{select(n),2};
            site_w = rf{select(n),3};
            site_h = rf{select(n),4};
            
            site_area = site_w * site_h;
            site_d{i,j}(n) = 2*sqrt(site_area/pi);
        end
        
        if isempty(select)
            mean_site_d{i,j} = [];
            subplot(length(xlocbin)-1,length(ylocbin)-1,i+(3-j)*(length(xlocbin)-1))  % j=3为第一排，j=2为第二排。。。
        else
            mean_site_d{i,j} = nanmean(site_d{i,j});
            subplot(length(xlocbin)-1,length(ylocbin)-1,i+(3-j)*(length(xlocbin)-1))  % j=3为第一排，j=2为第二排。。。
            site_hist{i,j} = hist(site_d{i,j},xbins);
            top = max(site_hist{i,j});
            hold on
            % 空值就不画了
            if sum(site_hist{i,j})~=0
                h_site{i,j} = bar(xbins,site_hist{i,j},1,'facecolor','w');
                plot(mean_site_d{i,j},top+1.5,'kv','markerfacecolor','k');
                text(mean_site_d{i,j}+9,top+1.5,num2str(mean_site_d{i,j}),'fontsize',10);
            else
                %axis off
            end
        end
        
        if i==2 && j==1
            xlabel('Equivalent diameter of receptive field (degree)');
        end
        if j==3 && i==1
            title('Grid x:16-17');
        elseif j==3 && i==2
            title('Grid x:18-19');
        elseif j==3 && i==3
            title('Grid x:20-21');
        end
        
        if i==1 && j==1
            ylabel('Grid y: 13-15');
        elseif i==1 && j==2
            ylabel('Grid y: 16-18');
        elseif i==1 && j==3
            ylabel('Grid y: 19-21');
        end
    end
end


%% ============ 2-D Visualization for big RF cell (Grid view from the top) =============== %
radius = 0.42;  % Radius of each hole (interval = 1)
% Plot grid outline
set(figure(103),'Position',[10 100 1100 700]); clf
hold on; 
axis equal ij;
gridRange = [10 23; 1 25]; % Grid plotting ragne = [xLow xHigh; yLow yHigh]
x1 = gridRange(1,1); x2 = gridRange(1,2);
y1 = gridRange(2,1); y2 = gridRange(2,2);
xlim([y1-2,y2+2]);

% Frame
interval = 5;
xLoc = intersect(x1:x2,0:interval:100);
yLoc = intersect(y1:y2,0:interval:100);
xLines = line(repmat([y1-1;y2+1],1,length(xLoc)),repmat(xLoc,2,1));  % xLines
yLines = line(repmat(yLoc,2,1), repmat([x1-1;x2+1],1,length(yLoc)));  % yLines
set(xLines,'LineWidth',5,'Color','g');
set(yLines,'LineWidth',5,'Color','g');
set(gca,'xtick', yLoc);
set(gca,'ytick', xLoc);

% Parallel drawing (to speed up)
xOffsets = repmat(x1:x2, y2-y1+1,1);
xOffsets = xOffsets(:)';
yOffsets = repmat([(y1:y2)+0.5*~mod(x1,2) (y1:y2)+0.5*mod(x1,2)], 1, fix((x2-x1 + 1)/2));
if mod(x2-x1+1,2) == 1
    yOffsets = [yOffsets y1:y2];
end
t = linspace(0,2*pi,100)';
xGrid = radius * repmat(sin(t),1,length(xOffsets)) + repmat(xOffsets,length(t),1);
yGrid = radius * repmat(cos(t),1,length(yOffsets)) + repmat(yOffsets,length(t),1);
set(fill(yGrid,xGrid,[1 1 1]),'LineWidth',1.5);    % White color

% Plot Big RF cell
c{1}='y';
c{2}='g';
c{3}='b';
c{4}='m';
c{5}='r';

for j = 5
    big_rf_cell = [rf{:,3}]>(j*10) & [rf{:,4}]>(j*10);
    big_rf = rf(big_rf_cell,:);
    big_rf_xloc = xloc(big_rf_cell);
    big_rf_yloc = yloc(big_rf_cell);
    
    for i= 1:length(big_rf)
        xCenter = big_rf_yloc(i);   %图上竖直方向为x，水平方向为y。。。
        yCenter = big_rf_xloc(i)+0.5 * ~mod(big_rf_yloc(i),2);
        % Paint
        t = linspace(0,2*pi,100)';
        xMap = radius * sin(t) + xCenter ;
        yMap = radius * cos(t) + yCenter;
        fill(yMap,xMap,c{j});
        alpha(0.5)
    end
end

%左脑
set(gca,'xdir','rev');


a=1;