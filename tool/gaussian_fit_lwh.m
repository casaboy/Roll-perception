% Lwh 2020.2
% Gaussian fit for gaussian distribution
function [a,b,c] = gaussian_fit_lwh(data_cum)
data_cum(isnan(data_cum(:,2)),:) = [];
data_cum = sort(data_cum);

% generate guess
q0 = ones(3,1);

% get a starting parameter estimate by testing a bunch of values
% y = a*exp(-((x-b)/c)^2) (不需要最后+base line)
% a range from [1 10 100 500]
% b range from [-10 0 10];
% c range from [0.01 0.1 1 10]
a_e = [1 10 100 500];
b_e = 0; % 默认为0
c_e = [0.01 0.1 1 10];
for i = 1:length(a_e)
%     for j = 1:length(b_e)
        for k = 1:length(c_e)
            q0(1,1) = a_e(i);
            q0(2,1) = b_e;
            q0(3,1) = c_e(k);
            errors(i,k) = cost_function(q0,data_cum);
        end
%     end
end
[min_indx1,min_indx3] = find(errors==min(min(errors)));
q0(1,1) = a_e(min_indx1(1));
q0(2,1) = b_e;
q0(3,1) = c_e(min_indx3(1));


% fo_ = fitoptions('method','NonlinearLeastSquares','Robust','On','Algorithm','T',...
%     'StartPoint',[q0(1,1) q0(2,1) q0(3,1)],'Lower',[0 -inf 0],'Upper',[Inf Inf Inf],'MaxFunEvals',5000);
% ft_ = fittype('a*exp(-((x-b)/c)^2)','dependent',{'y'},'independent',{'x'},'coefficients',{'a', 'b', 'c'});
% [cf_, goodness, output]= fit(data_cum(:,1),data_cum(:,2),ft_,fo_);

OPTIONS = optimset('MaxIter', 1e5,'MaxFunEvals', 1e5); % MaxIter允许迭代最大次数，MaxFunEvals：function估计最大数目（默认500）
% OPTIONS = optimset('PlotFcns',@optimplotfval);
[quick,fval,exitflag,output] = fminsearch(@(q)cost_function(q,data_cum),q0,OPTIONS);  % 算法：单纯性搜索法
% Output
a = quick(1,1);
b = quick(2,1);
c = quick(3,1);

function err = cost_function(q,data_cum)
TINY = 1e-10;
x = data_cum(:,1);
y = data_cum(:,2);

try
    n = data_cum(:,3);
catch
    n = ones(size(x));
end

z = q(1).*exp(-((x-q(2))./q(3)).^2);
z = z - (z > .999999)*TINY + (z < .0000001)*TINY;
% llik = n .* y .* log(z) +  n .* (1-y) .* log(1-z); % 交叉熵cross entropy, 一种常见的cost function
% err = -sum(llik);
err = norm(z-y);
