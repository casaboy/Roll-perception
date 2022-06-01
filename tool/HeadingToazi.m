function azi = HeadingToazi(x) % x in degree, azi in degree
azi = mod(90-x,360);  % Transferred from heading to azimuth