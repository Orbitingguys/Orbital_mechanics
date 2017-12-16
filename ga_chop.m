function D_v = ga_chop(x,orbital_parameters_1,orbital_parameters_2,mi)

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
mi = 398600;
[r_geo_1_0,v_geo_1_0] = kep2geo (orbital_parameters_1,mi);
[r_geo_2_0,v_geo_2_0] = kep2geo (orbital_parameters_2,mi);
options = odeset('Reltol',1e-13,'Abstol',1e-14);
X_1_0 = [r_geo_1_0;v_geo_1_0];
X_2_0 = [r_geo_2_0;v_geo_2_0];
[~,X_1] = ode113(@orbit_dynamics,linspace(0,x(1),3),X_1_0,options,mi);
[~,X_2] = ode113(@orbit_dynamics,linspace(0,x(2),3),X_2_0,options,mi);
tof = x(2)-x(1);

[~,~,~,~,v_transf_1,v_transf_2,t_p,~] = lambertMR(X_1(end,1:3),X_2(end,1:3),tof,mi,0,0,0,0);

if tof <= t_p
    D_v = NaN;
else
    D_v = norm(X_1(end,4:6)-v_transf_1)+norm(X_2(end,4:6)-v_transf_2);
end