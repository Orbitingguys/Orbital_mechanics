function [D_v,t_i,t_f,tof,k,New_min_D_v,i,j,vect_D_v,vect_diffD_v,X_1,X_2] = pork_chopWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi)

% TimeOption     .TypeMission:      -   if .TypeMission is 1 it means that you have a window of
%                                       launch so you need to define .t_i_max and also you have a
%                                       request for the arrival time which add a constrain that
%                                       has to be saved in matlab as .t_f_max
%
%                                   -   if .TypeMission is 2 you still have the launch window
%                                       but you don't deserve a constrain for the arrival time.
%
%                                   -   if .TypeMission is 3 you have to respect only an arrival
%                                       time
%
%                                   -   if .TypeMission is 4 you haven't any time constrain;
%   
%                .t_i_max:          -   Must be added for .TypeMission equal to 1 or 2           
%                                       for the other cases will be computed as the period of the
%                                       first orbit
%                                       
%                .t_f_max:          -   Must be added for .TypeMission equal to 1 or 3           
%                                       for the other cases will be computed as 4 times of the period 
%                                       of the second orbit 
%
%                WHEN A TIME IS NOT NEEDED JUST DON'T ADD IT AS A DATA!
%
%                .lenght_t_i        -   number of point examinated for the start orbit
%
%                .lenght_t_f        -   number of point examinated for the final orbit 
%
%                BE CAREFUL, WHIT THE Old TWO DATA, THE FUNCTION IS REALLY SLOW! 
%                

%% TIME SPAN

t_i_min = 0.01;                 % different from zero but close in order to prevent computing errors
t_f_min = 0.01;                 % It's impossible to reach the second orbit in 0 s but impossible cases will neglected by the code
%mi = 398600;
a_1 = orbital_parameters_1.a;
a_2 = orbital_parameters_2.a;
T_1 = 2*pi*sqrt(a_1^3/mi);    % period of the first orbit
T_2 = 2*pi*sqrt(a_2^3/mi);    % period of the second orbit


switch TimeOption.TypeMission
    
    
    case 1
        
        t_i_max = TimeOption.t_i_max;
        t_f_max = TimeOption.t_f_max;
        
    case 2
        
        t_i_max = TimeOption.t_i_max;
        t_f_max = 3*T_2;
        
    case 3
        
        t_i_max = T_1;
        t_f_max = TimeOption.t_f_max;
        
    case 4
        
       t_i_max = T_1;
       t_f_max = 3*T_2; 
       
end
    
 
%% ODE SETS

[r_geo_1_0,v_geo_1_0] = kep2geo (orbital_parameters_1,mi);
[r_geo_2_0,v_geo_2_0] = kep2geo (orbital_parameters_2,mi);
options = odeset('Reltol',1e-13,'Abstol',1e-14);
X_1_0 = [r_geo_1_0;v_geo_1_0];
X_2_0 = [r_geo_2_0;v_geo_2_0];


%% COMPUTING THE BEST DELTA_V

k = 1;
Old_min_D_v = 500;
diffD_v = 1;
stop = 0;
% m = TimeOption.lenght_t_i;
% n = TimeOption.lenght_t_f;
m = 100;
n = 100;
vect_D_v = [];
vect_diffD_v = [];
    
while  (k < 100   &&   diffD_v > 1e-3 && stop ~=1)
    k = k+1;
    [t_i,X_1] = ode113(@orbit_dynamics,linspace(t_i_min,t_i_max,m),X_1_0,options,mi);
    [t_f,X_2] = ode113(@orbit_dynamics,linspace(t_f_min,t_f_max,n),X_2_0,options,mi);
    D_v = zeros(m,n);
    tof = zeros(m,n);
    for i = 1:m
        for j = 1:n
            tof(i,j) = t_f(j)-t_i(i);
            if tof(i,j) <= 0
                tof(i,j) = NaN;
                D_v(i,j) = NaN;
            else
                [~,~,~,~,v_transf_1,v_transf_2,t_p,~] = lambertMR(X_1(i,1:3),X_2(j,1:3),tof(i,j),mi,0,0,0,0);
                if tof(i,j) <= t_p
                    tof(i,j) = NaN;
                    D_v(i,j) = NaN;
                else
                    D_v(i,j) = norm(X_1(i,4:6)-v_transf_1)+norm(X_2(j,4:6)-v_transf_2);
                end
            end
        end
    end
[~,i] = min(D_v);
[New_min_D_v,j] = min(min(D_v));
diffD_v = Old_min_D_v-New_min_D_v;
Old_min_D_v = New_min_D_v;
vect_diffD_v = [vect_diffD_v diffD_v];
vect_D_v = [vect_D_v New_min_D_v];

%% SEARCHING FOR STOP PARAMETER

if  i(j) ~= 1
    t_i_min = t_i(i(j)-1);
    X_1_0 = X_1(i(j)-1,:);
else
    t_i_min = t_i(i(j));
    X_1_0 = X_1(i(j),:);
end

if  i(j) ~= m
    t_i_max = t_i(i(j)+1);
else
    t_i_max = t_i(i(j));
end

if t_i_max == t_i_min
    stop = 1;
    fprintf('t_i_max = t_i_min\n')
end

if  j ~= 1 
    t_f_min = t_f(j-1);
    X_2_0 = X_2(j-1,:);
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
    fprintf('t_f_max = t_f_min\n')
end

if  abs(diffD_v) < 1e-3
    fprintf('Tollerance reached\n')
% elseif diffD_v < 0
%     error('merda');
%     m = m+50;
%     n = n+50;
%     k = 1;
%     Old_min_D_v = 500;
%     diffD_v = 1;
%     stop = 0;
end
    

end





