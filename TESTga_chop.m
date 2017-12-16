
clear;clc;close all
orbital_parameters_1.a = 12000;
orbital_parameters_1.ecc = 0;
orbital_parameters_1.theta = 120*pi/180;
orbital_parameters_1.RAAN = 0;
orbital_parameters_1.PA = 0;
orbital_parameters_1.INCLI = 0;
orbital_parameters_2.a = 9500;
orbital_parameters_2.ecc = 0.3;
orbital_parameters_2.theta = 250*pi/180;
orbital_parameters_2.RAAN = 0;
orbital_parameters_2.PA = 0;
orbital_parameters_2.INCLI = 0;
warning('off')
mi = 398600;
a_1 = orbital_parameters_1.a;
a_2 = orbital_parameters_2.a;
T_1 = 2*pi*sqrt(a_1^3/mi);    
T_2 = 2*pi*sqrt(a_2^3/mi);

[r_geo_1_0,v_geo_1_0] = kep2geo (orbital_parameters_1,mi);
[r_geo_2_0,v_geo_2_0] = kep2geo (orbital_parameters_2,mi);
options = odeset('Reltol',1e-13,'Abstol',1e-14);
X_1_0 = [r_geo_1_0;v_geo_1_0];
X_2_0 = [r_geo_2_0;v_geo_2_0];
% f = @(t_d,t_a) ga_chop(t_d,t_a,orbital_parameters_1,orbital_parameters_2,mi);

[x,fval] = ga(@ga_chop,2,[],[],[],[],[0,100000],[0,100000]);
% x = [5.8143e+03,1.5837e+04];
% D_v = ga_chop(x,orbital_parameters_1,orbital_parameters_2,mi);
