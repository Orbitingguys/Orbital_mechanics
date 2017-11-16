%% LAB 3
clear all;clc;close all
warning('off')
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
TimeOption.TypeMission = 4;
mi = 398600;
% TimeOption.lenght_t_i = 300;
% TimeOption.lenght_t_f = 300;

[D_v,t_i,t_f,tof,k,New_min_D_v,i,j,vect_D_v,vect_diff_D_v,X_1,X_2] = pork_chopWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi);

% surf(t_f,t_i,D_v,'EdgeColor','none');

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

TimeOption.TypeMission = 4;
[D_v,t_i,t_f,tof,k,New_min_D_v,i,j,vect_D_v,vect_diff_D_v,X_1,X_2] = pork_chopWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi);