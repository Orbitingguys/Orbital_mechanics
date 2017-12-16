



%% 1ST INTERPLANETARY LEG

clear; close all; clc

G = 6.67259e-20;
m_sun = 1.988919445342813e+030;
mi_sun = m_sun*G;

% DateSetup.body_1st = 2;
% DateSetup.body_2nd = 5;
% DateSetup.body_3rd = 6;
body_1st = 2;
body_2nd = 5;
body_3rd = 6;
% DateSetup.default = 0;
% DateSetup.TypeMission = 1;
date_d_min = [2018 2 31 0 0 0];
T0_1st = date2mjd2000(date_d_min)*3600*24;
mjd2000_0_1st = date2mjd2000(date_d_min);
t_a_min_1st = 0;
t_d_min_1st = 0;
orbital_parameters_venus_1st = OrbitalParameters(mjd2000_0_1st,body_1st);
a_venus = orbital_parameters_venus_1st.a;
orbital_parameters_jupiter_1st = OrbitalParameters(mjd2000_0_1st,body_2nd);
a_jupiter = orbital_parameters_jupiter_1st.a;
orbital_parameters_saturn_1st = OrbitalParameters(mjd2000_0_1st,body_3rd);
a_saturn = orbital_parameters_saturn_1st.a;
T_venus = 2*pi*sqrt(a_venus^3/mi_sun);
T_jupiter = 2*pi*sqrt(a_jupiter^3/mi_sun);
T_saturn = 2*pi*sqrt(a_saturn^3/mi_sun);
t_d_max_1st = T0_1st+T_venus;
% DateSetup.date_d_max = mjd20002date(t_d_max./(3600*24));
[r_geo_1_0,v_geo_1_0] = kep2geo (orbital_parameters_venus_1st,mi_sun);
[r_geo_2_0,v_geo_2_0] = kep2geo (orbital_parameters_jupiter_1st,mi_sun);
[r_geo_3_0,v_geo_3_0] = kep2geo (orbital_parameters_saturn_1st,mi_sun);
X_1_0_1st = [r_geo_1_0;v_geo_1_0];
X_2_0_1st = [r_geo_2_0;v_geo_2_0];
X_3_0_1st = [r_geo_3_0;v_geo_3_0];

mi_venus = astroConstants(12);
mi_jupiter = astroConstants(15);
mi_saturn = astroConstants(16);
% DateSetup.date_a_min = DateSetup.date_d_min;

t_a_max_1st = T0_1st+3*T_jupiter+T_venus;
% DateSetup.date_a_max = mjd20002date(t_a_max./(3600*24));

m = 10;
n = 30;
u = 60;


legtype_1st = 1;

[D_v_1st,t_d_1st,t_a_1st,tof_1st,X_1_1st,X_2_1st,v_transf_1st_1,v_transf_1st_2] = pork_chopLEG(X_1_0_1st,X_2_0_1st,...
    t_d_min_1st,t_d_max_1st,t_a_min_1st,t_a_max_1st,mi_sun,m,n,1);


%% 2ND INTERPLANETARY LEG

t_d_2nd = t_a_1st;
t_a_min_2nd = t_d_2nd(1);
t_a_max_2nd = t_d_2nd(1) + 3*T_saturn;
X_2_0_2nd = X_2_1st(1,:);
options = odeset('Reltol',1e-13,'Abstol',1e-14);
if t_d_2nd(1) ~= 0
    [~,X_3_2nd] = ode113(@orbit_dynamics,linspace(0,t_d_2nd(1),3),X_3_0_1st,options,mi_sun);
    X_3_0_2nd =  X_3_2nd(end,:);
else
    X_3_0_2nd = X_3_0_1st;
end


[D_v_2nd,t_d_2nd,t_a_2nd,tof_2nd,X_2_2nd,X_3_2nd,v_transf_2nd_1,v_transf_2nd_2] = pork_chopLEG(X_2_0_2nd,X_3_0_2nd,...
    t_d_2nd(1),t_d_2nd(end),t_a_min_2nd,t_a_max_2nd,mi_sun,n,u,2);

D_v_2nd = D_v_2nd';

[min_D_v_2,k] = min(D_v_2nd);

%% FLY BY

D_v_fby = zeros(m,u,n);
mjd2000_0_2nd = zeros(n,1);

for i = 1:n
    mjd2000_0_2nd(i) = t_d_2nd(i)/(60*60*24);
    orbital_parameters_jupiter_2nd = OrbitalParameters(mjd2000_0_2nd(i),body_2nd);
    for j = 1:m
        for w = 1:k
            v_1 = [v_transf_1st_2(j,i,1); v_transf_1st_2(j,i,2); v_transf_1st_2(j,i,3)];
            v_2 = [v_transf_2nd_1(i,w,1); v_transf_2nd_1(i,w,2); v_transf_2nd_1(i,w,3)];
            [~, ~, delta_v]  =  poweredflyby (v_1, v_2,orbital_parameters_jupiter_2nd,mi_jupiter);
            D_v_fby(j,w,i) = delta_v;
        end
    end
end

mD_v_fby = min(D_v_fby);
mmD_v_fby = min(mD_v_fby);
min_D_v_fby = min(mmD_v_fby)

% D_v_FB = zeros(n);
% for i = 2:n
%     [~,X_Flyby] = ode113(@orbit_dynamics,linspace(0,t_d_2nd(i),3),X_2_0_1st,options,mi_sun);
%     R_Flyby = X_Flyby(end,1:3); V_Flyby = X_Flyby(end,4:6);
%     orbital_parameters_planet = geo2kep(R_Flyby,V_Flyby,mi_sun);
%     for j = 1:m
%         v_1st_2 = [v_transf_1st_2(j,i,1);v_transf_1st_2(j,i,2);v_transf_1st_2(j,i,3)];
%         v_2nd_1 = [v_transf_2nd_1(i,k(i),1);v_transf_2nd_1(i,k(i),2);v_transf_2nd_1(i,k(i),3)];
%         if (norm(v_1st_2) && norm(v_2nd_1)) > 0
%             [~, ~, D_v_FB(i,j),deltaPA]  =  poweredflyby (v_1st_2,v_2nd_1,orbital_parameters_planet,mi_jupiter);
%         end
%         
%     end
% end

%% DEPARTURE FROM VENUS

%% ARRIVING AT SATURN

%% BEST DELTA V


% D_v = D_v_1st + min_D_v_2;