% 列联表卡方统计
% Lwh 20210324
% similar to function: crosstab
% add one-side test in fisher exact testing

function [p, chi2, df] = chi_square_Lwh(A,sided)

if nargin==1
    sided = 'yyds';
end

% "sided" only work in fisher exact testing!!!

A = round(A); % 确保都是整数

[r,c] = size(A);
sum_A = sum(sum(A));

% 计算理论频数T
T = zeros(r,c); % 观察频数
for i = 1:r
    for j = 1:c
        T(i,j) = sum(A(i,:))*sum(A(:,j))/sum_A;
    end
end

if  r== 2 && c == 2 % 四格表
    if sum_A >= 40 && sum(sum(T<5))==0 % 总样本量  n≥40，所有理论频数   T≥5 时，用Pearson卡方检验或似然比χ2 检验,他们结论基本一致。
        % Pearson卡方检验
        x_ = zeros(r,c); % 待求和的卡方分量
        for i = 1:r
            for j = 1:c
                x_(i,j)=(A(i,j)-T(i,j))^2 / T(i,j); % X2 = Σ(Oi-Ei)^2 / Ei
            end
        end
        chi2 = sum(sum(x_));
        df = (r-1)*(c-1);
        p = chi2pval(chi2,df);
        
    elseif sum_A >= 40 && sum(sum((T>=1 & T<5)))>0 % 总样本量  n≥40，但有理论频数1≤T＜5时，用连续性校正卡方（χ2）检验。
        % Yate's continuity correction
        x_ = zeros(r,c); % 待求和的卡方分量
        for i = 1:r
            for j = 1:c
                x_(i,j)=(abs(A(i,j)-T(i,j))-0.5)^2 / T(i,j); % X2 = Σ(|Oi-Ei| – 0.5)^2 / Ei
            end
        end
        chi2 = sum(sum(x_));
        df = (r-1)*(c-1);
        p = chi2pval(chi2,df);
        
    elseif sum_A < 40 || sum(sum(T<1))>0 % 总样本量  n＜40，或有理论频数  T＜1时，不能用卡方检验，用Fisher精确概率法检验(其实不属于卡方检验，因此没有卡方统计值)。
        % Fisher Exact testing
        if strcmpi(sided,'left')
            [~,p] = fishertest(A,'tail','left');
        elseif strcmpi(sided,'right')
            [~,p] = fishertest(A,'tail','right');
        else
            [~,p] = fishertest(A);
        end
        chi2 = nan;        
        % 当行数或者列数>2时，使用Fisher精确概率法检验的拓展:Freeman-Halton检验, 此时有卡方统计值
    else
        keybard
    end
else % R x C表
    if sum(sum(T<5))>(r*c/5) || sum(sum(T<1))>0
        keyboard
    else
        % Pearson卡方检验
        x_ = zeros(r,c); % 待求和的卡方分量
        for i = 1:r
            for j = 1:c
                x_(i,j)=(A(i,j)-T(i,j))^2 / T(i,j); % X2 = Σ(Oi-Ei)^2 / Ei
            end
        end
        chi2 = sum(sum(x_));
        df = (r-1)*(c-1);
        p = chi2pval(chi2,df);
    end
end

end