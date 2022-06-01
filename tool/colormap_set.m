% set color map
% example
% standard_color = [0 0 0; 1 0 0; 1 1 1];

function color_map_RGB = colormap_set(standard_color)

[color_num rgb_num] = size(standard_color);
if rgb_num~=3
    disp('Error Color Input!!!!!!!!!!!!')
    keyboard
end
if max(max(standard_color)) > 1
    disp('Error Color Input!!!!!!!!!!!!')
    keyboard
end

colorbin = 1000;
R_mat = [];
G_mat = [];
B_mat = [];

for i = 1:color_num-1
    R_mat = [R_mat, linspace(standard_color(i,1),standard_color(i+1,1),colorbin)];
    G_mat = [G_mat, linspace(standard_color(i,2),standard_color(i+1,2),colorbin)];
    B_mat = [B_mat, linspace(standard_color(i,3),standard_color(i+1,3),colorbin)];
end

color_map_RGB = [R_mat', G_mat', B_mat'];
end