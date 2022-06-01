function rot = axisToRotation(x,y,null_value)  % x: azimuth  y: rotation_axis
if (nargin < 3), null_value = -9999; end

for i = 1:length(x)
    if x(i) == 0 && y(i) == 0 % dot move from right to left
        rot(i) = 888;
    elseif x(i)==180 && y(i) ==0 % dot move from left to right
        rot(i) = 888;
    elseif x(i)==90 && y(i)==0 % Exp
        rot(i) = 0;
    elseif x(i)==270 && y(i)==0 % Con
        rot(i) = -180;
    elseif x(i)==0 && y(i)==90 % dot move CCW
        rot(i) = -90;
    elseif x(i)==0 && y(i)==270 % dot move CW
        rot(i) = 90;
    elseif x(i)==90 && y(i)==90 % dot move Exp+CCW
        rot(i) = -45;
    elseif x(i)==90 && y(i)==270 % dot move Exp+CW
        rot(i) = 45;
    elseif x(i)==270 && y(i)==90 % dot move Con+CCW
        rot(i) = -135;
    elseif x(i)==270 && y(i)==270 % dot move Con+CW
        rot(i) = 135;
    elseif x(i)== null_value
        rot(i) = null_value;
    else
        rot(i) = 888;
    end
end



