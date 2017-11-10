function [Local_min_D_V,R_d,R_a,V_d,V_a,v_transf_d,v_transf_a,ERROR] = WHILE_CHOP(t_i_min,t_i_max,t_f_min,t_f_max,m,n,mi,X_1_0,X_2_0,options)  

%% values for starting the while cycle

k = 1;                           % number of iteration
kmax = 100;                      % maximum number of iteration
toll = 1e-3;                     % tollerance accepted of the local minimum
Old_min_D_v = 500;               % value of the last minimum
diffD_v = 1;                     % difference between two consecutive minimum
stop = 0;                        % flag for the cycle
ERROR = zeros(kmax,1);           % pre-allocation of the memory
v_transf_1 = zeros(3,kmax,m,n);  % pre-allocation of the memory, do not remove, it's a 4D matrix!!
v_transf_2 = zeros(3,kmax,m,n);  % pre-allocation of the memory, do not remove, it's a 4D matrix!!

while  (k < kmax   &&   diffD_v > toll && stop ~=1)
        
k = k+1;   % the first iteration is dummy

[t_i,X_1] = ode113(@orbit_dynamics,linspace(t_i_min,t_i_max,m),X_1_0,options,mi); % ODE for the departure orbit with dynamics BCs
[t_f,X_2] = ode113(@orbit_dynamics,linspace(t_f_min,t_f_max,n),X_2_0,options,mi); % ODE for the departure orbit with dynamics BCs

D_v = zeros(m,n); % pre-allocation of the memory
tof = zeros(m,n); % pre-allocation of the memory

for i = 1:m                        % for every raw (t_i)
    for j = 1:n                    % for every column (t_f)
        tof(i,j) = t_f(j)-t_i(i);  % time of fligth
        if tof(i,j) <= 0
            tof(i,j) = NaN;        % putting NaN helps to plot the results
            D_v(i,j) = NaN;
        else
            [~,~,~,~,v_transf_1(:,k,i,j),v_transf_2(:,k,i,j),t_p,~] = lambertMR(X_1(i,1:3),X_2(j,1:3),tof(i,j),mi,0,0,0,0);
            if tof(i,j) <= t_p
                tof(i,j) = NaN;
                D_v(i,j) = NaN;
            else
                D_v(i,j) = norm(X_1(i,4:6)-(v_transf_1(:,k,i,j))')+norm(X_2(j,4:6)-(v_transf_2(:,k,i,j))'); % computing all Delta_v's combinations
            end
        end
    end
end

[~,i] = min(D_v);                   % Computing the vector of the raw position for every column minimum of D_v (as in the for cycle)
[New_min_D_v,j] = min(min(D_v));    % Saving the last local minimum
diffD_v = Old_min_D_v-New_min_D_v;  % Computing the difference between two consecutive minimum in order to reach the 'toll'
Old_min_D_v = New_min_D_v;          % Set to restart the while-cycle

%% SEARCHING FOR STOP PARAMETER
% As in the for cycle BCs are gonna be re-constrained

if  i(j) ~= 1
    t_i_min = t_i(i(j)-1);  % BC for the domine of the ODE of the departure orbit
    X_1_0 = X_1(i(j)-1,:);  % BC for the ODE of the departure orbit
else
    t_i_min = t_i(i(j));
    X_1_0 = X_1(i(j),:);
end



if  i(j) ~= m
    t_i_max = t_i(i(j)+1); % BC for the domine of the ODE of the departure orbit
else
    t_i_max = t_i(i(j));
end



if t_i_max == t_i_min
    stop = 1;
end



if  j ~= 1 
    t_f_min = t_f(j-1);
    X_2_0 = X_2(j-1,:);    % BC for the ODE of the arrival orbit
else
    t_f_min = t_f(j);
    X_2_0 = X_2(j,:);
end



if j ~= n
    t_f_max = t_f(j+1);
else
    t_f_max = t_f(j);
end



if t_f_max == t_f_min
    stop = 1;
end



if  abs(diffD_v) < 1e-3    % listing the ERROR matrix 
    ERROR(k) = 0;
elseif diffD_v < 0
    ERROR(k) = 1;
end

%% NEEDED VALUES

Local_min_D_V = New_min_D_v;           % Local computed minimum
v_transf_d = v_transf_1(:,k,i(j),j);   % Velocity vector in the transfer orbit for the local minimum conditions at the last iteration in the intersection whit the departure orbit
v_transf_a = v_transf_2(:,k,i(j),j);   % Velocity vector in the transfer orbit for the local minimum conditions at the last iteration in the intersection whit the arrival orbit
R_d = X_1(i(j),1:3);                   % Radius of the first intersection for the local minimum conditions
V_d = X_1(i(j),4:6);                   % Velocity vector in the departure orbit for the local minimum conditions at the last iteration
R_a = X_2(j,1:3);                      % Radius of the second intersection for the local minimum conditions
V_a = X_2(j,4:6);                      % Velocity vector in the arrival orbit for the local minimum conditions at the last iteration

end
