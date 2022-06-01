% Multi-unit spike rate adjustment
% 消除随着时间的变化：如随着时间FR变大 or 变小  (空值或0：默认不作处理；1:detrend处理, 2多项式拟合，3 local regression)
% Lwh
% e.g.
% spike_rates = spike_rate_adjust(spike_rates,1);

function spike_rates = spike_rate_adjust(spike_rates_raw,methodCode)

if nargin==1
    methodCode=0;
end

if methodCode~=0
    switch methodCode
        case 1
            % 去掉随时间的趋势(去掉线性趋势)
            residual = detrend(spike_rates_raw);
            trend = spike_rates_raw - residual;
            
            mu = mean(spike_rates_raw);
            spike_rates = residual + mu;
            
%             clf
%             figure(1)
%             subplot(1,2,1)
%             hold on
%             plot(spike_rates_raw);
%             plot(trend,':r');
%             plot(residual,'m');
%             plot(zeros(size(residual)),':k')
%             legend('Original Data','Trend','Detrended Data',...
%                    'Mean of Detrended Data','Location','northwest')
%             ylimt = ylim;
%             
%             subplot(1,2,2)
%             plot(spike_rates)
%             ylim(ylimt);

        case 2
            disp('Not complete');
            keyboard
            %         % 去掉非线性趋势 (1-3次多项式)
            %         m=3;
            %         [p s] = polyfit([1:length(spike_rates_raw)],spike_rates_raw,m);
            %         x1 = [1:0.1:length(spike_rates_raw)];
            %         y1 = polyval(p,x1);
            
        case 3
            % local regression
            x = 1:length(spike_rates_raw);
            y = spike_rates_raw;
            r = ksrlin(x,y);  % 修改N值为trial数目
            residual = spike_rates_raw - r.f;
            mu = mean(spike_rates_raw);
            spike_rates = residual + mu;   
    end
else
    spike_rates = spike_rates_raw;
end
end