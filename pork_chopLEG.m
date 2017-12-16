function [D_v,t_d,t_a,tof,X_1,X_2,v_transf_1,v_transf_2] = pork_chopLEG(X_1_0,X_2_0,t_d_min,t_d_max,t_a_min,t_a_max,mi,length_t_d,length_t_a,legtype)



options = odeset('Reltol',1e-13,'Abstol',1e-14);

[t_d,X_1] = ode113(@orbit_dynamics,linspace(t_d_min,t_d_max,length_t_d),X_1_0,options,mi);
[t_a,X_2] = ode113(@orbit_dynamics,linspace(t_a_min,t_a_max,length_t_a),X_2_0,options,mi);

v_transf_1 = zeros(length_t_d,length_t_a,3);
v_transf_2 = zeros(length_t_d,length_t_a,3);
D_v = zeros(length_t_d,length_t_a);
tof = zeros(length_t_d,length_t_a);
for i = 1:length_t_d

    for j = 1:length_t_a
        tof(i,j) = t_a(j)-t_d(i);
        if tof(i,j) <= 0
            tof(i,j) = NaN;
            D_v(i,j) = NaN;
        else 
           [~,~,~,~,v_transf_1(i,j,:),v_transf_2(i,j,:),t_p,~] = lambertMR(X_1(i,1:3),X_2(j,1:3),tof(i,j),mi,0,0,0,0);
           v_1 = [v_transf_1(i,j,1),v_transf_1(i,j,2),v_transf_1(i,j,3)];
           v_2 = [v_transf_2(i,j,1),v_transf_2(i,j,2),v_transf_2(i,j,3)];
           if tof(i,j) <= t_p
               tof(i,j) = NaN;
               D_v(i,j) = NaN;
           else  
               if legtype == 1
                    D_v(i,j) = norm(X_1(i,4:6)-v_1);
                elseif legtype == 2
                    D_v(i,j) = norm(X_2(j,4:6)-v_2);
                elseif legtype == 0
                    D_v(i,j) = norm(X_1(i,4:6)-v_1)+norm(X_2(j,4:6)-v_2);
                end
           end
        end
    end
end




