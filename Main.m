options = odeset('Reltol',1e-13,'Abstol',1e-14);
R_I = [-5.5106881000000003e+03  1.5889450000000002e+03  5.8383016000000007e+03];
V_I = [-3.8639999999999999e+00 -6.5410000000000004e+00 -1.2840000000000000e+00];
t_0 = 0.1;
t_f = 100000;
mi = 398600;
% orbital_parameter.a = ;
% orbital_parameter.ecc = ;
% orbital_parameter.PA = ;
% orbital_parameter.RAAN = ;
% orbital_parameter.INCLI = ;
% orbital_parameter.theta = ;
equat_Rt = 6378.1363;
polar_Rt = 6356.7523142;
X_0 = [R_I(1); R_I(2); R_I(3); V_I(1); V_I(2); V_I(3)];

I = imread('earth.jpg'); 
[x_t, y_t, z_t] = ellipsoid (0,0,0,equat_Rt, equat_Rt, polar_Rt);
terra = surf(x_t, y_t, -z_t,'Edgecolor', 'none');
set(terra,'FaceColor','texturemap','Cdata',I)
% set(gca,'Color','none','visible','off')
view(127.5,30)
axis equal


[T,X] = ode113(@orbit_dynamics,[t_0 t_f],X_0,options,mi);
hold on
plot3(X(:,1),X(:,2),X(:,3))

