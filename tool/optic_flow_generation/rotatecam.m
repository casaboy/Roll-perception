function [x_,y_,z_] = rotatecam(x,y,z,rotate_ang)

y_ = y;
x_ = x*cos(rotate_ang)-z*sin(rotate_ang);
z_ = z*cos(rotate_ang)+x*sin(rotate_ang);