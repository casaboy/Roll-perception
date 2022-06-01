% poisson spike generator (firing rate固定，不随着时间变化)
% Lwh 20200506
% ti+1 = ti - ln(xrand)/r
% ti: spike time (interspike interval)
% xrand = uniformly in the range between 0 and 1
% r; firing rate


t0 = 0;
% 进行100个trial
for m = 1:100
    xrand = rand(1,100);
    r = 10; % 设定firing rate
    for i = 1:50 % 前50的spike
        if i == 1
            t(m,i) = t0-log(xrand(i))/r;
        else
            t(m,i) = t(m,i-1)-log(xrand(i-1))/r;
        end
    end
    
    r = 100; % 设定firing rate
    for i = 51:100 % 后50个spike
        t(m,i) = t(m,i-1)-log(xrand(i-1))/r;
    end
end

spike = ones(1,100);
figure
plot(t(1,1:50),spike(1,1:50),'k.');
hold on
plot(t(1,51:100),spike(1,51:100),'r.');

% 每段都是50个spike，所以firing rate低的时间更长，firing rate高的时间更短
% 计算实际给出的firing rate：r = n/T
for m = 1:100
    r_es1(m) = 50/(t(m,50)-t(m,1));
    r_es2(m) = 50/(t(m,100)-t(m,51));
end

mean_spike1 = mean(r_es1);
mean_spike2 = mean(r_es2);
variance_spike1 = std(r_es1)^2; % 不太准确，理论上mean=variance
variance_spike2 = std(r_es2)^2;

fano_factor1 = variance_spike1/mean_spike1;
fano_factor2 = variance_spike2/mean_spike2;


aaa =1;