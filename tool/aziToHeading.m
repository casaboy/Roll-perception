function azi = aziToHeading(x,null_value)
if (nargin < 2), null_value = -9999; end

azi = nan(size(x));
for i = 1:length(x)
    if x(i) == 888
        azi(i) = 888; % add by Lwh 20201103
    elseif x(i) == null_value
        azi(i) = null_value; % add by Lwh 20201218
    else
        azi(i) = mod(270-x(i),360)-180; % Transferred from azimuth to heading
    end
end
