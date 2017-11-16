function [t]=timesinceper(theta,h,ecc,mu)
%This function returns the time passed from the periapse till the possition
%of true anomally theta on a given orbit in its plane
f=@(x)(h^3./(mu.^2*(1+ecc*cos(x)).^2));
t=quadgk(f,0,theta,'RelTol',1e-10,'AbsTol',1e-14);