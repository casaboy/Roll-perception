function [col_str] = num2ExcelName(num_loc) % Convert xls column number to column name (such as "A", "AB", ...)
test = 2;
old = 0;
x = 0;
while test >= 1
    old = 26^x + old;
    test = num_loc/old;
    x = x + 1;
end
num_letters = x - 1;
str_array = zeros(1,num_letters);
for i = 1:num_letters
    loc = floor(num_loc/(26^(num_letters-i)));
    num_loc = num_loc - (loc*26^(num_letters-i));
    str_array(i) = char(65 + (loc - 1));
end
col_str = strcat(str_array(1:length(str_array)));
end