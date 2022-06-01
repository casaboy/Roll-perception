% angle between two 2 dim vector, Lwh 202010
% x1: two dimension vector or column x*2
% y1: two dimension vector
function y1 = cossAngle(x1,x2)
x1 = squeeze(x1);
x2 = squeeze(x2);

if size(x1,2) == 2 && size(x2,2) == 2 && size(x1,1) ~=1 && size(x2,1) == 1
    for i = 1:size(x1,1)
        x = x1(i,:);
        y = dot(x,x2)/(norm(x)*norm(x2));
        y1(i,1) = acosd(y); % in to degree
    end
else
    try
        y = dot(x1,x2)/(norm(x1)*norm(x2));
        y1 = acosd(y); % in to degree
    catch
        disp('Error Input!!!!!!!!!!!!')
        keyboard
    end
end
