% Lwh
% used for calculate P vaule for microstimulation through probit regression
% respond to extreme behavior in animals to obtain an acurate P value 
% e.g.
% 1. monkey do not make any choice under a specific condition
%    for example, in heaidng -6 degree, monkey choose all CW or CCW, then
%    for the translation task, this consition will be NaN (Left+Right = 0)
% 2. monkey choose only one side in the microstimulation
%    for example, monkey choose all right target in the translation task.
% 3. 除了极端的1至2个condition，全部选一边

function [X,Y] = AdjustXY_for_probit(x,yctrl,ystim, repctrl, repstim)
% x: unique_degree
% yctrl: Psy_correct, percentage of corret non-stimulated trial
% ystim: Psy_correct, percentage of corret stimulated trial
% repctrl: repetition of non-stimulated trial 
% repstim: repetition of timulated trial 

% check data
if size(x,2)~=1
    x = x';
end

if size(yctrl,2)~=1
    yctrl = yctrl';
end

if size(ystim,2)~=1
    ystim = ystim';
end

if size(repctrl,2)~=1
    repctrl = repctrl';
end

if size(repstim,2)~=1
    repstim = repstim';
end


%% 1. remove the NaN condition
if sum(isnan(ystim))>0 || sum(isnan(yctrl))>0
    no_choice = any([isnan(ystim),isnan(yctrl)],2);
    % remove this condition
    yctrl(no_choice) = [];
    ystim(no_choice) = [];
    x(no_choice) = [];
    repctrl(no_choice) = [];
    repstim(no_choice) = [];
end

%% 2. all choose one side: 手动在原有unique_degree外足够远的位置增加一个假设的选择，使其选择为1或0
if sum(ystim == 1) == length(x) % 在左边增加一个选择0
    ystim = [0; ystim];
    yctrl = [0; yctrl];
    x = [min(x)*2; x];
elseif sum(ystim == 0) == length(x) % 在右边增加一个选择1
    ystim = [ystim; 1];
    yctrl = [yctrl; 1];
    x = [x;max(x)*2];
end

%% 3. 除了极端的1至2个condition，全部选一边: 100%的个数>总degree数-2，或者0%的个数大于总degree数-2
% 通过插值增加数据点（WARNING：不知道对真实值的估计有没有影响。。。。不过一般这个时候的p值都非常非常显著，所以p值的绝对数值应该对结果是没影响的）
if sum(ystim == 1) >= length(x)-2 || sum(ystim == 0) >= length(x)-2
    xx = [min(x):0.2:max(x)]'; % 暂时都用0.2插值
    yyctrl = interp1(x,yctrl,xx,'pchip'); % spline interpolate 不合适
    yystim = interp1(x,ystim,xx,'pchip');
    
    yyctrl(yyctrl<0) = 0; % 插值得到的Y可能会比0要小，调整为0
    yystim(yystim<0) = 0;
    
    yyctrl(yyctrl>1) = 1; % 插值得到的Y可能会比1要大，调整为1
    yystim(yystim>1) = 1;
    
    % 此时repetition全部假定为repctrl、repstim的均值
    rep_temp = mean([repctrl;repstim]);
    repctrl = ones(size(yyctrl)).* rep_temp;
    repstim = ones(size(yystim)).* rep_temp;
else
    xx = x;
    yyctrl = yctrl;
    yystim = ystim;
end

% data for probit regression
X = [];
X(:,1) = [xx;xx];
X(:,2) = [zeros(size(xx)); ones(size(xx))]; % 0: not stim; 1: stim
X(:,3) = X(:,1).* X(:,2);

Y = [];
Y(:,1) = [yyctrl.* repctrl; yystim .* repstim]; % correct trial
Y(:,2) = [repctrl; repstim]; % correct trial
