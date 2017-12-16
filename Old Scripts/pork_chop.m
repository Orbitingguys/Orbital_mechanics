function [D_v,t_i,t_f,tof,min_D_v,i,j,X_1,X_2,t_p] = pork_chop(orbital_parameters_1,orbital_parameters_2,t_i_min,t_i_max,t_f_min,t_f_max,mi,length_t_i,length_t_f)

tic

[r_geo_1_0,v_geo_1_0] = kep2geo (orbital_parameters_1,mi);
[r_geo_2_0,v_geo_2_0] = kep2geo (orbital_parameters_2,mi);
options = odeset('Reltol',1e-13,'Abstol',1e-14);
X_1_0 = [r_geo_1_0;v_geo_1_0];
X_2_0 = [r_geo_2_0;v_geo_2_0];
[t_i,X_1] = ode113(@orbit_dynamics,linspace(t_i_min,t_i_max,length_t_i),X_1_0,options,mi);
[t_f,X_2] = ode113(@orbit_dynamics,linspace(t_f_min,t_f_max,length_t_f),X_2_0,options,mi);
v_transf_1 = zeros(3,length_t_i,length_t_f);
D_v = zeros(length_t_i,length_t_f);
tof = zeros(length_t_i,length_t_f);
t_p = zeros(length_t_i,length_t_f);
for i = 1:length_t_i
    for j = 1:length_t_f
        tof(i,j) = t_f(j)-t_i(i);
        if tof(i,j) <= 0
            tof(i,j) = NaN;
            D_v(i,j) = NaN;
        else 
           [~,~,~,~,v_transf_1(:,i,j),v_transf_2,t_p(i,j),~] = lambertMR(X_1(i,1:3),X_2(j,1:3),tof(i,j),mi,0,0,0,0);
           if tof(i,j) <= t_p(i,j)
               tof(i,j) = NaN;
               D_v(i,j) = NaN;
           else
               D_v(i,j) = norm(X_1(i,4:6)-(v_transf_1(:,i,j))')+norm(X_2(j,4:6)-v_transf_2);
           end
        end
    end
end

figure
surf(t_f,t_i,D_v,'EdgeColor','none');
figure
contour(t_f,t_i,D_v,length_t_i);

[~,i] = min(D_v);
[min_D_v,j] = min(min(D_v));

[orbital_parameters_transf] = geo2kep(X_1(i(j),1:3),(v_transf_1(:,i(j),j))',mi);
a_transf = orbital_parameters_transf.a;
T_transf = 2*pi*sqrt(a_transf^3/mi); 
t_i_min = 0.001;
t_i_max = T_transf;
X_transf_0 = [X_1(i(j),1:3),(v_transf_1(:,i(j),j))'];
[~,X_transf] = ode113(@orbit_dynamics,[t_i_min t_i_max],X_transf_0,options,mi);

t_i_max = 2*pi*sqrt(orbital_parameters_1.a^3/mi);
t_f_max = 2*pi*sqrt(orbital_parameters_2.a^3/mi);

[~,X_1] = ode113(@orbit_dynamics,linspace(t_i_min,t_i_max,length_t_i),X_1_0,options,mi);
[~,X_2] = ode113(@orbit_dynamics,linspace(t_f_min,t_f_max,length_t_f),X_2_0,options,mi);

figure
hold on
plot3(X_transf(:,1),X_transf(:,2),X_transf(:,3))
plot3(X_1(:,1),X_1(:,2),X_1(:,3))
plot3(X_2(:,1),X_2(:,2),X_2(:,3))
hold off

toc
        
