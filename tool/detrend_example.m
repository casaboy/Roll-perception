sig = [0 1 -2 1 0 1 -2 1 0];      % signal with no linear trend
trend = [0 1 2 3 4 3 2 1 0];      % two-segment linear trend
x = sig+trend;                    % signal with added trend
y = detrend(x,'linear',5)         % breakpoint at 5th element

figure(3)
plot(sig,'k-');
hold on
plot(trend,'b-');
