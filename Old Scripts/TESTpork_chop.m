%% LAB 3

clear all;clc;close all
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

[D_v,t_i,t_f,tof,min_D_v,i,j,X_1,X_2] = pork_chop(orbital_parameters_1,orbital_parameters_2,0.01,T_1,0.01,3*T_2,mi,200,200);




%% MARS ORBIT

clear all;clc;close all
warning('off')
mi=1.32712440018E11 ;
% t_d=linspace(0,200*24*3600,100);
% t_a=linspace(150*24*3600,400*24*3600,100);
orbital_parameters_1.ecc=1.730508176494217E-02;
orbital_parameters_1.a=1.497209137830794E+08;
orbital_parameters_1.RAAN=degtorad(1.613676165270718E+02);
orbital_parameters_1.PA=degtorad(3.034472453387018E+02);
orbital_parameters_1.INCLI=degtorad(1.141631465551398E-03);
orbital_parameters_1.theta=degtorad(5.513047086433101E+01);

orbital_parameters_2.ecc=9.354121433212964E-02;
orbital_parameters_2.a=2.279359554264199E+08;
orbital_parameters_2.RAAN=degtorad(4.954660277763613E+01);
orbital_parameters_2.PA=degtorad(2.865236216522881E+02);
orbital_parameters_2.INCLI=degtorad(1.849354955586735E+00);
orbital_parameters_2.theta=degtorad(2.538807159694742E+02);

a_1 = orbital_parameters_1.a;
a_2 = orbital_parameters_2.a;
T_1 = 2*pi*sqrt(a_1^3/mi);    
T_2 = 2*pi*sqrt(a_2^3/mi);



[D_v,t_i,t_f,tof,min_D_v,i,j,X_1,X_2,t_p] = pork_chop(orbital_parameters_1,orbital_parameters_2,0,200*24*3600,150*24*3600,400*24*3600,mi,200,200);
% [D_v,t_i,t_f,tof,min_D_v,i,j,X_1,X_2] = pork_chop(orbital_parameters_1,orbital_parameters_2,0.01,T_1,0.01,3*T_2,mi,100,100);
