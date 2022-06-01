function [theta,psi] = camproject(x,y,z)
% x_ = x-campos(1);
% y_ = y-campos(2);
% z_ = z-campos(3);

[theta,psi,r] = cart2sph(x,y,z);
i0 = r<10;
i1 = theta<pi&theta>0;
psi(psi>pi) = -(2*pi-psi(psi>pi));
% i2 = (psi<pi/4&psi>0)|(psi>3/2*pi&psi<2*pi);
i2 = (psi>-pi/2)&(psi<pi/2);
i = i0&i1&i2;
theta = theta(i);
psi = psi(i);

