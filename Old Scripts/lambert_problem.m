function [v_i,v_f] = lambert_problem(r_i,r_f,delta_t) 

% P_i and P_f are identified by r_i and r_f
% This script computes the time of fligth between P_i and P_f;
norm_r_i = norm(r_i);
norm_r_f = norm(r_f);
delta_theta = asin(dot(r_i,r_f)/(norm_r_i*norm_r_f));
c = sqrt(norm_r_i^2 + norm_r_f^2 - 2*norm_r_i*norm_r_f*cos(delta_theta)); % chord between P_i and P_f
p = (norm_r_i+norm_r_f+c)/2;     %semiperimetre
 
a_m = p/2;
alpha = 2*asin(sqrt(p/(2*a_m)));
if delta_theta > 0 && delta_theta < pi
    beta = 2*asin(sqrt((p-c)/(2*a_m)));
else 
    beta = pi-2*asin(sqrt((p-c)/(2*a_m)));
end


delta_t_m = (alpha-beta-sin(alpha)+sin(beta))/sqrt(mi/a_m^3);

