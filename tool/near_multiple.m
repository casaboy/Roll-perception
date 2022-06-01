% 取一个数最接近的3/4/5的倍数
% 取一个数最接近的3/4的倍数
% Lwh  20210904 for xbin in ploting
function m = near_multiple(n,multi)

time_multuple = [];
for i = 1:length(multi)
    time_multuple(i) = ceil(n/multi(i))*multi(i);
end

% find the nearst multiple and bigger than n
m = min(time_multuple);

if m < n
    keyboard
end

end