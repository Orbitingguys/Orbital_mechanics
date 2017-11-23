%% LAB 5
clc;clear all;close all
v_s_minus = [31.5; 4.69; 0];
v_s_plus = [38.58; 0; 0];
R_vect = [0; -150000000; 0];
R = norm(R_vect);
G = 6.67259e-20;
msun = 1.988919445342813e+030;
mi_sun = msun*G;
mi_earth = 398600;
v_circ = sqrt(mi_sun/R);
v_circ_vect = [v_circ; 0; 0 ];
v_inf_minus_vect = v_s_minus - v_circ_vect;
gamma = atan(v_inf_minus_vect(2)/v_inf_minus_vect(1));
v_inf_minus = norm(v_inf_minus_vect);
a_inc = - mi_earth/v_inf_minus^2;
v_inf_plus_vect = v_s_plus - v_circ_vect;
v_inf_plus = norm(v_inf_plus_vect);
a_out = - mi_earth/v_inf_plus^2;
delta = gamma;

x_0 = [1.1 1.1 pi/3 pi/3 1000]';

x = fsolve(@(x)hyp_system(x,v_inf_plus,v_inf_minus,delta,mi_earth),x_0);

e_inc = x(1);
e_out = x(2);
delta_inc = x(3);
delta_out = x(4);
r_p = x(5);

v_p_inc = sqrt(mi_earth/(a_inc*(1-e_inc^2)))*(1+e_inc);
v_p_out = sqrt(mi_earth/(a_out*(1-e_out^2)))*(1+e_out);
delta_v = v_p_out-v_p_inc;
theta_inf_inc = acos(-1/e_inc);

if sin(theta_inf_inc) ~= sqrt((e_inc^2 - 1)/e_inc)
    theta_inf_inc = 2*pi - theta_inf_inc;
end

orbital_parameters.a = a_inc;
orbital_parameters.ecc = e_inc;
orbital_parameters.RAAN = 0;
orbital_parameters.PA = 0;
orbital_parameters.INCLI = 0;
orbital_parameters.theta_inf = theta_inf_inc;
orbital_parameters.incoming = true;
[~,r_inc_vect,~] = hyp_orbit(orbital_parameters,mi_earth);

% theta_inc = linspace(-theta_inf_inc,0,10000);
% r_inc = zeros(length(theta_inc),1);
% r_inc_vect = zeros(length(theta_inc),2);
% 
% for i = 1:length(theta_inc)
%     
% r_inc(i) = a_inc*(1-e_inc^2)./(1+e_inc.*cos(theta_inc(i)));
% r_inc_vect(i,:) = [r_inc(i)*cos(theta_inc(i)) r_inc(i)*sin(theta_inc(i))];
% 
% end
plot(r_inc_vect(1,8000:end),r_inc_vect(2,8000:end),'k')
hold on
theta_inf_out = acos(-1/e_out);

if sin(theta_inf_out) ~= sqrt((e_out^2 - 1)/e_out)
    theta_inf_out = 2*pi - theta_inf_out;
end


orbital_parameters.a = a_out;
orbital_parameters.ecc = e_out;
orbital_parameters.RAAN = 0;
orbital_parameters.PA = 0;
orbital_parameters.INCLI = 0;
orbital_parameters.theta_inf = theta_inf_out;
orbital_parameters.incoming = false;
[~,r_out_vect,~] = hyp_orbit(orbital_parameters,mi_earth);


% theta_out = linspace(0,theta_inf_out,10000);
% r_out = zeros(length(theta_out),1);
% r_out_vect = zeros(length(theta_out),2);
% 
% for i = 1:length(theta_out)
%     
% r_out(i) = a_out*(1-e_out^2)/(1+e_out*cos(theta_out(i)));
% r_out_vect(i,:) = [r_out(i)*cos(theta_out(i)) r_out(i)*sin(theta_out(i))];
% 
% end
plot(r_out_vect(1,1:2000),r_out_vect(2,1:2000),'y')



