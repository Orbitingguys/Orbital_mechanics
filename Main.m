clc; clear all; close all
options = odeset('Reltol',1e-13,'Abstol',1e-14);
R_I = [-5.51e+03  1.58e+03  5.83e+03];
V_I = [-3.86e+00 -6.54e+00 -1.28e+00];
t_0 = 0.01;
% t_f = 100000;
mi = 398600;
equat_Rt = 6378.1363;
polar_Rt = 6356.7523142;
X_0 = [R_I(1); R_I(2); R_I(3); V_I(1); V_I(2); V_I(3)];

% I = imread('earth.jpg'); 
% [x_t, y_t, z_t] = ellipsoid (0,0,0,equat_Rt, equat_Rt, polar_Rt);
% terra = surf(x_t, y_t, -z_t,'Edgecolor', 'none');
% set(terra,'FaceColor','texturemap','Cdata',I)
% % set(gca,'Color','none','visible','off')
% view(127.5,30)
% axis equal

[orbital_parameters] = geo2kep(R_I,V_I,mi);
T = 2*pi*sqrt(orbital_parameters.a^3/mi);
t_f = 4*T;
[t,X] = ode113(@orbit_dynamics,[t_0 t_f],X_0,options,mi);
V_geo = [X(:,4),X(:,5),X(:,6)];
R_geo = [X(:,1),X(:,2),X(:,3)];
% hold on
% plot3(X(:,1),X(:,2),X(:,3))

omega_planet = 2*pi/(24*3600);
[map_tif] = imread('NE1_LR_LC_SR_W_DR.tif');
map_tfw = worldfileread('NE1_LR_LC_SR_W_DR.tfw');
figure
mapshow(map_tif, map_tfw);
hold on
[alpha,delta,lambda,phi,index] = groundtrack(R_geo,omega_planet,t);
n_index = length(index);
plot(lambda(1:index(1)),phi(1:index(1)),'k')

for i = 2:n_index
plot(lambda((index(i-1)+1):index(i)),phi((index(i-1)+1):index(i)),'k')
end
plot(lambda(index(end)+1:end),phi(index(end)+1:end),'k');
plot(lambda(1),phi(1),'b.','MarkerSize',20)
plot(lambda(end),phi(end),'r.','MarkerSize',20)

hold off

