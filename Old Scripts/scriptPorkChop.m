clear all;clc;close all
warning('off')
mi=1.32712440018E11 ;
t_d=linspace(0,200*24*3600,50);
t_a=linspace(0,400*24*3600,50);
orbital_parameters_1.ecc=1.730508176494217E-02;
orbital_parameters_1.a=1.497209137830794E+08;
orbital_parameters_1.RAAN=degtorad(1.613676165270718E+02);
orbital_parameters_1.PA=degtorad(3.034472453387018E+02);
orbital_parameters_1.INCLI=degtorad(1.141631465551398E-03);
orbital_parameters_1.theta=degtorad(5.513047086433101E+01);
[R01,V01] = kep2geo (orbital_parameters_1,mi,orbital_parameters_1.theta);
orbital_parameters_2.ecc=9.354121433212964E-02;
orbital_parameters_2.a=2.279359554264199E+08;
orbital_parameters_2.RAAN=degtorad(4.954660277763613E+01);
orbital_parameters_2.PA=degtorad(2.865236216522881E+02);
orbital_parameters_2.INCLI=degtorad(1.849354955586735E+00);
orbital_parameters_2.theta=degtorad(2.538807159694742E+02);
[R02,V02] = kep2geo (orbital_parameters_2,mi,orbital_parameters_2.theta);

[Dv,ERROR,TPAR,TOF]=porkchop(t_d,t_a,mi,orbital_parameters_1,orbital_parameters_2);