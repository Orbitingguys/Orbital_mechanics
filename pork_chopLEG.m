function [D_v,t_i,t_f,tof,X_1,X_2,v_transf] = pork_chopLEG(orbital_parameters_1,orbital_parameters_2,t_i_min,t_i_max,t_f_min,t_f_max,mi,length_t_i,length_t_f,legtype)


[r_geo_1_0,v_geo_1_0] = kep2geo (orbital_parameters_1,mi);
[r_geo_2_0,v_geo_2_0] = kep2geo (orbital_parameters_2,mi);
options = odeset('Reltol',1e-13,'Abstol',1e-14);
X_1_0 = [r_geo_1_0;v_geo_1_0];
X_2_0 = [r_geo_2_0;v_geo_2_0];
[t_i,X_1] = ode113(@orbit_dynamics,linspace(t_i_min,t_i_max,length_t_i),X_1_0,options,mi);
[t_f,X_2] = ode113(@orbit_dynamics,linspace(t_f_min,t_f_max,length_t_f),X_2_0,options,mi);
v_transf = zeros(3,length_t_i,length_t_f);
v_transf_1 = zeros(3,length_t_i,length_t_f);
v_transf_2 = zeros(3,length_t_i,length_t_f);
D_v = zeros(length_t_i,length_t_f);
tof = zeros(length_t_i,length_t_f);
for i = 1:length_t_i
    for j = 1:length_t_f
        tof(i,j) = t_f(j)-t_i(i);
        if tof(i,j) <= 0
            tof(i,j) = NaN;
            D_v(i,j) = NaN;
        else 
           [~,~,~,~,v_transf_1(:,i,j),v_transf_2(:,i,j),t_p,~] = lambertMR(X_1(i,1:3),X_2(j,1:3),tof(i,j),mi,0,0,0,0);
           if tof(i,j) <= t_p
               tof(i,j) = NaN;
               D_v(i,j) = NaN;
           else
               if legtype == 1
                    D_v(i,j) = norm(X_1(i,4:6)-v_transf_1(:,i,j)');
                    v_transf(:,i,j) = v_transf_1(:,i,j);
                elseif legtype == 2
                    D_v(i,j) = norm(X_2(j,4:6)-v_transf_2(:,i,j)');
                    v_transf(:,i,j) = v_transf_2(:,i,j);
                elseif legtype == 0
                    D_v(i,j) = norm(X_1(i,4:6)-v_transf_1(:,i,j)')+norm(X_2(j,4:6)-v_transf_2(:,i,j)');
                    v_transf = 0;
                end
           end
        end
    end
end




