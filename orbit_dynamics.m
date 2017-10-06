%% Lab 01 

function [dX] = orbit_dynamics(~,X,mi)


r_x = X(1);
r_y = X(2);
r_z = X(3);
v_x = X(4);
v_y = X(5);
v_z = X(6);

dr_x = v_x;
dr_y = v_y;
dr_z = v_z;
R = [r_x,r_y,r_z];
norm_R = norm(R);
dv_x = -mi/norm_R^3 * r_x;
dv_y = -mi/norm_R^3 * r_y;
dv_z = -mi/norm_R^3 * r_z;

dX=[dr_x, dr_y, dr_z, dv_x, dv_y, dv_z]';

end
