function [r,ppp,h] = plot_corr_line(xx,yy,varargin)
% Lwh 202004
% Lwh 202101
% Examples
% [r_square,ppp,h] = plot_corr_line(xx,yy,'MethodOfCorr','pearson','FittingMethod',2,'LineStyles','r-','LineWidth',2);
% plot_corr_line(maxFR_group{i}(cell_type==ct,1),maxFR_group{i}(cell_type==ct,2),'MethodOfCorr','pearson','FittingMethod',2,'LineStyles',{'-','color',c{ct}},'LineWidth',2);

% Validate the type parameter, and set a handle to the correlation function.
% ------ Parse input parameters -------
paras = inputParser;
addOptional(paras,'MethodOfCorr','Pearson');
addOptional(paras,'FittingMethod',1); % Type 1: Normal squared error (最小二乘法LSM，竖直距离最小); Type 2: Perpendicular error (solved by PCA，垂直距离最小)
addOptional(paras,'LineStyles','k-'); % 每一条拟合线的linestyle
addOptional(paras,'LineWidth',1); % 每一条拟合线的linestyle
parse(paras,varargin{:});

lineProps = paras.Results.LineStyles;
if ~iscell(lineProps), lineProps={lineProps}; end
colorThis = find(strcmp('color',lineProps));

% xx,yy中任意一项有NAN则去除
noxx = isnan(xx);
noyy = isnan(yy);
remove_xy = logical(noxx | noyy);
xx = xx(~remove_xy);
yy = yy(~remove_xy);

% remove Inf/-Inf if have
noxx2 = logical(xx == inf | xx == -inf);
noyy2 = logical(yy == inf | yy == -inf);
remove_xy2 = logical(noxx2 | noyy2);
xx = xx(~remove_xy2);
yy = yy(~remove_xy2);


if ~isempty(xx)
    [rrr,ppp] = corr(xx,yy,'type',paras.Results.MethodOfCorr);
    r = rrr;
else
    r = nan;
    ppp = nan;
end

% if isnan(ppp)
%     disp('********** Something WRONG in corr **********');
%     keyboard;
% end

if ~isempty(xx)
    if paras.Results.FittingMethod == 2 % Minimize perpendicular distance (PCA)
        % New version of regressing perpendicular error HH20180613
        fitType1 = regress_perp(xx,yy,1,0.05);  % regress_perp(x,y,option,alpha)   option: 1--free intercept (default)  2--set intercept 0  3--set slope 2 (???)
        if ~isempty(fitType1)
            linPara(1) = fitType1.k;
            linPara(2) = fitType1.b;
        else
            linPara = [nan nan];
        end
        linParaSE = [std(fitType1.bootK) std(fitType1.bootB)];
        %     h.fitType1 = fitType1;
        
        % -- Plotting
        xxx = linspace(min(xx),max(xx),150);
        Y = linPara(1) * xxx + linPara(2);
        xxx = xxx(min(yy) <= Y & Y <= max(yy));
        Y = Y(min(yy) <= Y & Y <= max(yy));
        
    elseif paras.Results.FittingMethod == 1 % LMS method
        [b,~,stat]= glmfit(xx,yy); % glmfit can directly give us std of coeffs (se)
        linPara = fliplr(b'); % Note the definitions different
        linParaSE = fliplr(stat.se');
        
        % -- Plotting
        xxx = linspace(min(xx),max(xx),100);
        Y = polyval(linPara,xxx);
    end
    
    % -- Saving
    hold on
    if ppp<0.05
        if length(lineProps)==1 % e.g. 'r:'
            h = plot(xxx,Y,lineProps{1},'linewidth',paras.Results.LineWidth);
        else % e.g. '-','color',c{1}
            if ~isempty(colorThis)
                h = plot(xxx,Y,lineProps{1},'color',lineProps{colorThis+1},'linewidth',paras.Results.LineWidth);
            else
                h = plot(xxx,Y,lineProps{1},'color','k','linewidth',paras.Results.LineWidth);
            end
        end
    else % p>0.05 写死了线的格式
        if ~isempty(colorThis)
            h = plot(xxx,Y,'color',lineProps{colorThis+1},'linestyle','-','linewidth',paras.Results.LineWidth);
        else
            h = plot(xxx,Y,'color','k','linestyle','-','linewidth',paras.Results.LineWidth);
        end
    end
end
end