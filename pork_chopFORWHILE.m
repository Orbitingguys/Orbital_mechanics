function [Local_min_D_v,Total_min_D_v,ERROR] = pork_chopFORWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi)



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
%                                       for the other cases will be computed as 3 times of the period 
%                                       of the second orbit 
%
%                WHEN A TIME IS NOT NEEDED JUST DON'T ADD IT AS A DATA!                

tic

%% TIME SETS

t_i_min = 0.001;                 % different from zero but close in order to prevent computing errors
t_f_min = 0.001;                 % It's impossible to reach the second orbit in the very same time as the departure but impossible cases will neglected by the script
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

%% FIRST PORK CHOP CYCLE

    m = 200;      % number of steps in which the departure time will be discretized
    n = 200;      % number of steps in which the arrival time will be discretized
    
    % dynamic equation integration
    
    [t_i,X_1] = ode113(@orbit_dynamics,linspace(t_i_min,t_i_max,m),X_1_0,options,mi); % integration along the setted time for the departure orbit 
    [t_f,X_2] = ode113(@orbit_dynamics,linspace(t_f_min,t_f_max,n),X_2_0,options,mi); % integration along the setted time for the arrival orbit 
    [~,X_1_first] = ode113(@orbit_dynamics,linspace(t_i_min,T_1,m),X_1_0,options,mi); % integration along the whole time for the departure orbit
    [~,X_2_first] = ode113(@orbit_dynamics,linspace(t_f_min,T_2,n),X_2_0,options,mi); % integration along the whole time for the arrival orbit
     
    
    D_v = zeros(m,n); % pre-allocation of the memory
    tof = zeros(m,n); % pre-allocation of the memory
    
    for i = 1:m                       % for each raw (t_i)
        for j = 1:n                   % for each column (t_f)
            tof(i,j) = t_f(j)-t_i(i); % time of fligth
            if tof(i,j) <= 0    
                tof(i,j) = NaN;       % putting NaN helps to plot the results
                D_v(i,j) = NaN;
            else
                [~,~,~,~,v_transf_1,v_transf_2,t_p,~] = lambertMR(X_1(i,1:3),X_2(j,1:3),tof(i,j),mi,0,0,0,0);
                if tof(i,j) <= t_p
                    tof(i,j) = NaN;
                    D_v(i,j) = NaN;
                else
                    D_v(i,j) = norm(X_1(i,4:6)-v_transf_1)+norm(X_2(j,4:6)-v_transf_2); % computing all Delta_v's combinations
                end
            end
        end
    end
    
%% PlOTS

% 3D surface plot

figure                                      
surf(t_f,t_i,D_v,'EdgeColor','none');
title('Pork Chop Surface','FontSize',13)
xlabel('arrival time')
ylabel('departure time')
zlabel('Delta_v')

% contour plot

figure                                      
contour(t_f,t_i,D_v,m);
title('Pork Chop Contour','FontSize',13)
xlabel('arrival time')
ylabel('departure time')


%% COMPUTING THE BEST DELTA_V
% the aim of the for cycle is to run a while pork chop cycle inside around some local
% minimum points (found at the first pork chop)in order to compute the best global delta_v

% pre-allocation of the memory

Local_min_D_v = zeros(n,1);
ERROR = zeros(100,n);
R_d = zeros(3,n);
R_a = zeros(3,n);
V_d = zeros(3,n);
V_a = zeros(3,n);
v_transf_d = zeros(3,n);
v_transf_a = zeros(3,n);

[~,h] = min(D_v); % [minimum,index] = min(matrix) gives us a vector (minimum) which has in the k-th position the minimum of the ... 
                  % ... k-th matrix column and index is the vector in which presents the raw of the matrix for every minimum 

                  
% in order to introduce the while cycle we need to see how the BCs change inside the for cycle

for z = 1:n       % for cycle for every element of the minimum vector of the matrix (the last one computed)

    if h(z) > 3
        t_i_min = t_i(h(z)-3);   % BC for the domine of the ODE of the departure orbit
        X_1_0 = X_1(h(z)-3,:);   % BC for the ODE of the departure orbit
    else
        t_i_min = t_i(1);
        X_1_0 = X_1(1,:);        
    end
    
    if h(z) < m-3
        t_i_max = t_i(h(z)+3);   % BC for the domine of the ODE of the departure orbit
    else
        t_i_max = t_i(end);
    end
    
    if z > 2
        t_f_min = t_f(z-2);      % BC for the domine of the ODE of the arrival orbit
        X_2_0 = X_2(z-2,:);      % BC for the ODE of the arrival orbit
    else
        t_f_min = t_f(1);
        X_2_0 = X_2(1,:);
    end
    
    if z < n-2
        t_f_max = t_f(z+2);      % BC for the domine of the ODE of the arrival orbit
    else
        t_f_max = t_f(end);
    end
    
    [Local_min_D_v(z),R_d(:,z),R_a(:,z),V_d(:,z),V_a(:,z),v_transf_d(:,z),v_transf_a(:,z),err] = WHILE_CHOP(t_i_min,t_i_max,t_f_min,t_f_max,round(m/8),round(n/8),mi,X_1_0,X_2_0,options);
    
    ERROR(:,z) = err; % assembling the error matrix
    
    
   
end

[Total_min_D_v,f] = min(Local_min_D_v); % Computing the global minimum!

%% PRINTING INTERESTING VALUES

fprintf('\n\n')
fprintf('BEST DELTA_V: %g [km/s] \n\n\n',Total_min_D_v)
fprintf('DEPARTURE RADIUS:\n\n')
fprintf('x: %g [km] \n',R_d(1,f))
fprintf('y: %g [km] \n',R_d(2,f))
fprintf('z: %g [km] \n\n\n',R_d(3,f))
fprintf('ARRIVAL RADIUS:\n\n')
fprintf('x: %g [km] \n',R_a(1,f))
fprintf('y: %g [km] \n',R_a(2,f))
fprintf('z: %g [km] \n\n\n',R_a(3,f))
fprintf('DEPARTURE VELOCITY @departure orbit: \n\n')
fprintf('x: %g [km/s] \n',V_d(1,f))
fprintf('y: %g [km/s] \n',V_d(2,f))
fprintf('z: %g [km/s] \n\n\n',V_d(3,f))
fprintf('DEPARTURE VELOCITY @transfer orbit: \n\n')
fprintf('x: %g [km/s] \n',v_transf_d(1,f))
fprintf('y: %g [km/s] \n',v_transf_d(2,f))
fprintf('z: %g [km/s] \n\n\n',v_transf_d(3,f))
fprintf('ARRIVAL VELOCITY @arrival orbit:\n\n')
fprintf('x: %g [km/s] \n',V_a(1,f))
fprintf('y: %g [km/s] \n',V_a(2,f))
fprintf('z: %g [km/s] \n\n\n',V_a(3,f))
fprintf('ARRIVAL VELOCITY @transfer orbit:\n\n')
fprintf('x: %g [km/s] \n',v_transf_a(1,f))
fprintf('y: %g [km/s] \n',v_transf_a(2,f))
fprintf('z: %g [km/s] \n',v_transf_a(3,f))


%% ORBIT REPRESENTATION (geo coordinate plot)
[orbital_parameters_transf] = geo2kep(R_d(:,f),v_transf_d(:,f),mi);
a_transf = orbital_parameters_transf.a;
T_transf = 2*pi*sqrt(a_transf^3/mi); 
t_i_min = 0.001;
t_i_max = T_transf;
X_transf_0 = [R_d(:,f),v_transf_d(:,f)];
[~,X_transf] = ode113(@orbit_dynamics,[t_i_min t_i_max],X_transf_0,options,mi);

figure
hold on
start = plot3(X_transf(1,1),X_transf(1,2),X_transf(1,3),'o','MarkerEdgeColor','none','MarkerFaceColor','green','MarkerSize',8);
stop = plot3(X_2_first(end,1),X_2_first(end,2),X_2_first(end,3),'o','MarkerEdgeColor','none','MarkerFaceColor','red','MarkerSize',8);
trasf = plot3(X_transf(:,1),X_transf(:,2),X_transf(:,3),'y');
departure = plot3(X_1_first(:,1),X_1_first(:,2),X_1_first(:,3),'g');
arrival = plot3(X_2_first(:,1),X_2_first(:,2),X_2_first(:,3),'r');
xlabel('x')
ylabel('y')
zlabel('z')
title('Orbit Representation','FontSize',13)
legend([trasf,departure,arrival,start,stop],{'transfer orbit','departure orbit','arrival orbit','start position for the departure orbit', ...
    'start position for the arrival orbit'})


hold off


toc  
