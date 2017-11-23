

%% DEPARTURE FROM VENUS

%% 1st INTERPLANETARY LEG

clear all; close all; clc
warning('off');
G = 6.67259e-20;
msun = 1.988919445342813e+030;
mi = msun*G;

DateSetup.ibody = 2;
DateSetup.jbody = 5;
DateSetup.default = 0;
DateSetup.TypeMission = 1;

orbital_parameters_venus= OrbitalParameters(1,DateSetup.ibody);
a_venus = orbital_parameters_venus.a;
orbital_parameters_jupiter = OrbitalParameters(1,DateSetup.jbody);
a_jupiter = orbital_parameters_jupiter.a;

T_venus = 2*pi*sqrt(a_venus^3/mi);
T_jupiter = 2*pi*sqrt(a_jupiter^3/mi);

DateSetup.date_d_min = [2017 2 31 0 0 0];
T0 = date2mjd2000(DateSetup.date_d_min)*3600*24;

T_d_max = T0+T_venus;
DateSetup.date_d_max = mjd20002date(T_d_max./(3600*24));


DateSetup.date_a_min = DateSetup.date_d_min;

T_a_max = T0+2*T_jupiter+T_venus;
DateSetup.date_a_max = mjd20002date(T_a_max./(3600*24));
[printed_value_leg1,ERROR_leg1] = PORKCHOP_PROCEDURE(DateSetup);

return

%% 2ND INTERPLANETARY LEG

DateSetup.ibody = 5;
DateSetup.jbody = 6;
DateSetup.default = 0;
[printed_value_leg2,ERROR_leg2] = PORKCHOP_PROCEDURE(DateSetup);

%% ARRIVING AT jupiter