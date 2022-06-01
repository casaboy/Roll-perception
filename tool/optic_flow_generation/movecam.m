function [x_,y_,z_] = movecam(x,y,z,move_direction,move_distance)

v_ = move_direction/norm(move_direction);
x_ = x-v_(1)*move_distance;
y_ = y-v_(2)*move_distance;
z_ = z-v_(3)*move_distance;