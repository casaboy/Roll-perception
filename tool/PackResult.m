function result = PackResult(varargin)
% Pack variables into stucture result with their respective names
% HH20141124
for i = 1:length(varargin)
    try
    result.(inputname(i)) = varargin{i}; % Dynamic structure fields
    catch
        keyboard
    end
end