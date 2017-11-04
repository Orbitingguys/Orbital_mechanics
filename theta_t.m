function [theta]=theta_t(t,mu,h,ecc)
%This function returns the true anomaly value for a body on a given orbit
%in his plane, giving as an imput the time since periapsis (time from the
%moment the body pass through periapsis)
a=h^2/(mu*(1-ecc^2));
T=2*pi*sqrt(a^3/mu);
Me=2*pi*t/T; %We calculate the mean anomaly Me
E=fsolve(@(E)(E-ecc*sin(E)-Me),0); %We calculate the eccentric anomaly
theta=2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2)); %We get the true anomaly, but this value will be included in the interval [-pi,pi] 
if theta<0
    theta=2*pi+theta; %We define the true anomaly in the right interval, [0,2pi]
end
