function [printed_value,ERROR] = pork_chopPARPOLFORWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi,T0,m,n)

% TimeOption     .TypeMission:      -   if .TypeMission is 1 it means that you have a window of
%                                       launch so you need to define .t_i_max and also you have a
%                                       request for the arrival time which add a constrain that
%                                       has to be saved in matlab as .t_f_min,.t_f_max
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
%                .t_f_min, t_f_max: -   Must be added for .TypeMission equal to 1 or 3           
%                                       for the other cases will be computed as 3 times of the period 
%                                       of the second orbit 
%
%                The initial time at the departure orbit is automatically
%                computed as zero
%
%                WHEN A TIME IS NOT NEEDED JUST DON'T ADD IT AS A DATA!                

tic

%% TIME SETS

t_i_min = 0.001;                 % different from zero but close in order to prevent computing errors
a_1 = orbital_parameters_1.a;
a_2 = orbital_parameters_2.a;
T_1 = 2*pi*sqrt(a_1^3/mi);       % period of the first orbit
T_2 = 2*pi*sqrt(a_2^3/mi);       % period of the second orbit


switch TimeOption.TypeMission 
    
    
    case 1
        
        t_i_max = TimeOption.t_i_max;
        t_f_min = TimeOption.t_f_min;
        t_f_max = TimeOption.t_f_max;
        
    case 2

        t_i_max = TimeOption.t_i_max;
        t_f_min = 0.001;                 % It's impossible to reach the second orbit in the very same time as the departure but impossible cases will neglected by the script
        t_f_max = 3*T_2;
        
    case 3

        t_i_max = T_1;
        t_f_min = TimeOption.t_f_min;
        t_f_max = TimeOption.t_f_max;
        
    case 4
       
        t_i_max = T_1;
        t_f_min = 0.001;                 % It's impossible to reach the second orbit in the very same time as the departure but impossible cases will neglected by the script
        t_f_max = 3*T_2; 
       
end
    
 
%% ODE SETS

% m = 200;                                                    % number of steps in which the departure time will be discretized
% n = 200;                                                    % number of steps in which the arrival time will be discretized

[r_geo_1_0,v_geo_1_0] = kep2geo (orbital_parameters_1,mi);
[r_geo_2_0,v_geo_2_0] = kep2geo (orbital_parameters_2,mi);
options = odeset('Reltol',1e-13,'Abstol',1e-14);            % option set for the ODE integration

X_1_0 = [r_geo_1_0;v_geo_1_0];                              % starting condition for the ODE integration @ departure orbit (at t=0 [s])
X_2_0 = [r_geo_2_0;v_geo_2_0];                              % starting condition for the ODE integration @ arrival orbit (at t=0 [s])
X_0_departure = X_1_0;                                      % in order to be used later
X_0_arrival = X_2_0;                                        % in order to be used later

[~,X_departure_orbit] = ode113(@orbit_dynamics,linspace(t_i_min,T_1,m),X_1_0,options,mi); % integration along the whole time for the departure orbit
[~,X_arrival_orbit] = ode113(@orbit_dynamics,linspace(0.001,T_2,n),X_2_0,options,mi);     % integration along the whole time for the arrival orbit


switch TimeOption.TypeMission
    
    case 1
        
        if t_f_min ~= 0
        [~,X_2] = ode113(@orbit_dynamics,[0 t_f_min],X_2_0,options,mi); % integration needed in order to obtain the radius of the condition related to t_f_min in the arrival orbit
        X_2_0 = X_2(end,:);                                             % radius and velocity at the time t_f_min @ the arrival orbit
        end
    case 2
        
                                                                        % start conditions are fine as saved before   
        
    case 3
        
        if t_f_min ~= 0
        [~,X_2] = ode113(@orbit_dynamics,[0 t_f_min],X_2_0,options,mi); % integration needed in order to obtain the radius of the condition related to t_f_min in the arrival orbit
        X_2_0 = X_2(end,:);                                             % radius and velocity at the time t_f_min @ the arrival orbit
        end
    case 4
        
                                                                        % start conditions are fine as saved before
  
end


%% FIRST PORK CHOP CYCLE

    
    % dynamic equation integration
    
    [t_i,X_1] = ode113(@orbit_dynamics,linspace(t_i_min,t_i_max,m),X_1_0,options,mi); % integration along the setted time for the departure orbit 
    [t_f,X_2] = ode113(@orbit_dynamics,linspace(t_f_min,t_f_max,n),X_2_0,options,mi); % integration along the setted time for the arrival orbit 
     
    
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

%In order to represent the dates in the axes we need to get their tau, for that we first get the
%mjd2000 vector of times:

mjd2000_d = (t_i+T0)./(3600*24);
mjd2000_a = (t_f+T0)./(3600*24);

%Now we can define the tau for each time of both vectors
tau_d = zeros(length(mjd2000_d),1);
tau_a = zeros(length(mjd2000_a),1);

for i = 1:length(mjd2000_d)
    tau_d(i) = mjd20002tau(mjd2000_d(i));
end
for j=1:length(mjd2000_a)
    tau_a(j) = mjd20002tau(mjd2000_a(j));
end

% 3D SURFACE PLOT

figure                                      
surf(tau_d,tau_a,D_v','EdgeColor','none');

%Settings for surface plot
title('Pork Chop Surface','FontSize',13)
xlabel('Arrival time')
ylabel('Departure time')
zlabel('Delta_v')
tickDeparture = linspace(tau_d(1),tau_d(end),20);
tickArrival = linspace(tau_a(1),tau_a(end),20);
set(gca,'xtick',tickDeparture,'XTickLabelRotation',45);
set(gca,'ytick',tickArrival,'YTickLabelRotation',45);
datetick('x',1,'keepticks');
datetick('y',1,'keepticks');

% CONTOUR PLOT

[vt] = valuesTOF(mjd2000_d,mjd2000_a);
mjd2000of = tof./(3600*24);
figure                                      
contour(tau_d,tau_a,D_v',50);
hold on
contour(tau_d,tau_a,mjd2000of',vt,'LineColor','k','ShowText','on');

%Settings for contour plot
colorbar;
caxis([0,100]);
title('Pork Chop Contour','FontSize',13)
xlabel('Arrival time')
ylabel('Departure time')
set(gca,'xtick',tickDeparture,'XTickLabelRotation',45);
set(gca,'ytick',tickArrival,'YTickLabelRotation',45);
datetick('x',1,'keepticks');
datetick('y',1,'keepticks');
colorbar;
hold off

%% COMPUTING THE BEST DELTA_V
% the aim of the for cycle is to run a while pork chop cycle inside around some local
% minimum points (found at the first pork chop)in order to compute the best global delta_v

% pre-allocation of the memory

ERROR = zeros(100,n);
Local_min_D_v = zeros(n,1);
R_maneouvre_1 = zeros(3,n);
R_maneouvre_2 = zeros(3,n);
V_maneouvre_1 = zeros(3,n);
V_maneouvre_2 = zeros(3,n);
v_transf_man_1 = zeros(3,n);
v_transf_man_2 = zeros(3,n);
t_maneouvre_1 = zeros(1,n); 
t_maneouvre_2 = zeros(1,n); 

[~,h] = min(D_v); % [minimum,index] = min(matrix) gives us a vector (minimum) which has in the k-th position the minimum of the ... 
                  % ... k-th matrix column and index is the vector in which presents the raw of the matrix for every minimum 

                  
% in order to introduce the while cycle we need to see how the BCs change inside the for cycle

parfor_progress(n);
parpool;
parfor z = 1:n       % for cycle for every element of the minimum vector of the matrix (the last one computed)

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
   
    [Local_min_D_v(z),R_maneouvre_1(:,z),R_maneouvre_2(:,z),V_maneouvre_1(:,z),V_maneouvre_2(:,z),v_transf_man_1(:,z),v_transf_man_2(:,z),t_maneouvre_1(z),t_maneouvre_2(z),err] = WHILE_CHOP(t_i_min,t_i_max,t_f_min,t_f_max,round(m/20),round(n/20),mi,X_1_0,X_2_0,options);
    
    ERROR(:,z) = err;            % assembling the error matrix
    
    parfor_progress;
end

delete(gcp('nocreate'));
[Total_min_D_v,f] = min(Local_min_D_v); % Computing the global minimum! the minimum position is the f-th

[~,X_dep_t02Tman2] = ode113(@orbit_dynamics,linspace(0,t_maneouvre_2(f),10),X_0_departure,options,mi); % in order to obtain the position of the planet 1 @maneouvre time 2
[~,X_arr_t02Tman1] = ode113(@orbit_dynamics,linspace(0,t_maneouvre_1(f),10),X_0_arrival,options,mi);   % in order to obtain the position of the planet 2 @maneouvre time 1

% in order to transform the times of the maneuvers in date format ([yyyy mm dd hh mm ss]):

mjd2000_man1 = (t_maneouvre_1(f)+T0)./(3600*24);
mjd2000_man2 = (t_maneouvre_2(f)+T0)./(3600*24);
days_of_flight = mjd2000_man2-mjd2000_man1;

Date_man1 = mjd20002date(mjd2000_man1);
Date_man2 = mjd20002date(mjd2000_man2);

Datestr_man1 = datestr(Date_man1);
Datestr_man2 = datestr(Date_man2);

%% PRINTING INTERESTING VALUES

fprintf('\n\n')
fprintf('BEST DELTA_V: %g [km/s] \n\n\n',Total_min_D_v)
fprintf('DEPARTURE RADIUS @departure orbit:\n\n')
fprintf('x: %g [km] \n',X_departure_orbit(1,1))
fprintf('y: %g [km] \n',X_departure_orbit(1,2))
fprintf('z: %g [km] \n\n\n',X_departure_orbit(1,3))
fprintf('ARRIVAL RADIUS @arrival orbit:\n\n')
fprintf('x: %g [km] \n',X_arrival_orbit(1,1))
fprintf('y: %g [km] \n',X_arrival_orbit(1,2))
fprintf('z: %g [km] \n\n\n',X_arrival_orbit(1,3))
fprintf('DEPARTURE RADIUS @transfer orbit(first maneouvre point):\n\n')
fprintf('x: %g [km] \n',R_maneouvre_1(1,f))
fprintf('y: %g [km] \n',R_maneouvre_1(2,f))
fprintf('z: %g [km] \n\n\n',R_maneouvre_1(3,f))
fprintf('ARRIVAL RADIUS @transfer orbit(second maneouvre point):\n\n')
fprintf('x: %g [km] \n',R_maneouvre_2(1,f))
fprintf('y: %g [km] \n',R_maneouvre_2(2,f))
fprintf('z: %g [km] \n\n\n',R_maneouvre_2(3,f))
fprintf('DEPARTURE RADIUS @second maneouvre point:\n\n')
fprintf('x: %g [km] \n',X_dep_t02Tman2(end,1))
fprintf('y: %g [km] \n',X_dep_t02Tman2(end,2))
fprintf('z: %g [km] \n\n\n',X_dep_t02Tman2(end,3))
fprintf('ARRIVAL RADIUS @first maneouvre point:\n\n')
fprintf('x: %g [km] \n',X_arr_t02Tman1(end,1))
fprintf('y: %g [km] \n',X_arr_t02Tman1(end,2))
fprintf('z: %g [km] \n\n\n',X_arr_t02Tman1(end,3))
fprintf('DEPARTURE VELOCITY @departure orbit: \n\n')
fprintf('x: %g [km/s] \n',V_maneouvre_1(1,f))
fprintf('y: %g [km/s] \n',V_maneouvre_1(2,f))
fprintf('z: %g [km/s] \n\n\n',V_maneouvre_1(3,f))
fprintf('DEPARTURE VELOCITY @transfer orbit: \n\n')
fprintf('x: %g [km/s] \n',v_transf_man_1(1,f))
fprintf('y: %g [km/s] \n',v_transf_man_1(2,f))
fprintf('z: %g [km/s] \n\n\n',v_transf_man_1(3,f))
fprintf('ARRIVAL VELOCITY @arrival orbit:\n\n')
fprintf('x: %g [km/s] \n',V_maneouvre_2(1,f))
fprintf('y: %g [km/s] \n',V_maneouvre_2(2,f))
fprintf('z: %g [km/s] \n\n\n',V_maneouvre_2(3,f))
fprintf('ARRIVAL VELOCITY @transfer orbit:\n\n')
fprintf('x: %g [km/s] \n',v_transf_man_2(1,f))
fprintf('y: %g [km/s] \n',v_transf_man_2(2,f))
fprintf('z: %g [km/s] \n\n\n',v_transf_man_2(3,f))
fprintf('TIME @first maneouvre point: \n\n')
fprintf('t: %g [s] \n\n\n',t_maneouvre_1(f))
fprintf('TIME @second maneouvre point: \n\n')
fprintf('t: %g [s] \n\n\n',t_maneouvre_2(f))
fprintf('DATE @first maneouvre point: \n\n')
fprintf('%s \n\n\n',Datestr_man1)
fprintf('DATE @second maneouvre point: \n\n')
fprintf('%s \n\n\n',Datestr_man2)
fprintf('DAYS OF FLIGHT for the Express Mission: \n\n')
fprintf('~%g [days] \n\n\n',round(days_of_flight))


printed_value.best_Dv = Total_min_D_v;
printed_value.StartRad_DepOrbit = [X_departure_orbit(1,1),X_departure_orbit(1,2),X_departure_orbit(1,3)];
printed_value.StartRad_ArrOrbit = [X_arrival_orbit(1,1),X_arrival_orbit(1,2),X_arrival_orbit(1,3)];
printed_value.RadMan_1 = [R_maneouvre_1(1,f),R_maneouvre_1(2,f),R_maneouvre_1(3,f)];
printed_value.RadMan_2 = [R_maneouvre_2(1,f),R_maneouvre_2(2,f),R_maneouvre_2(3,f)];
printed_value.Man2Rad_DepOrbit = [X_dep_t02Tman2(end,1),X_dep_t02Tman2(end,3),X_dep_t02Tman2(end,3)];
printed_value.Man1Rad_ArrOrbit = [X_arr_t02Tman1(end,1),X_arr_t02Tman1(end,2),X_arr_t02Tman1(end,3)];
printed_value.VelDepOrbit_Man_1 = [V_maneouvre_1(1,f),V_maneouvre_1(2,f),V_maneouvre_1(3,f)];
printed_value.VelTransfOrbit_Man_1 = [v_transf_man_1(1,f),v_transf_man_1(3,f),v_transf_man_1(3,f)];
printed_value.VelArrOrbit_Man_2 = [V_maneouvre_2(1,f),V_maneouvre_2(2,f),V_maneouvre_2(3,f)];
printed_value.VelTransfOrbit_Man_2 = [v_transf_man_2(1,f),v_transf_man_2(3,f),v_transf_man_2(3,f)];
printed_value.Time_Man_1 = t_maneouvre_1(f);
printed_value.Time_Man_2 = t_maneouvre_2(f);
printed_value.Date_Man_1 = Date_man1;
printed_value.Date_Man_2 = Date_man2;

%% ORBIT REPRESENTATION (geo coordinate plot)

[orbital_parameters_transf] = geo2kep(R_maneouvre_1(:,f),v_transf_man_1(:,f),mi);
a_transf = orbital_parameters_transf.a;
T_transf = 2*pi*sqrt(a_transf^3/mi); 
t_i_min = 0.001;
t_i_max = T_transf;
X_transf_0 = [R_maneouvre_1(:,f),v_transf_man_1(:,f)];
[~,X_transf] = ode113(@orbit_dynamics,[t_i_min t_i_max],X_transf_0,options,mi);



figure
hold on

StartRad_DepOrbit = plot3(X_departure_orbit(1,1),X_departure_orbit(1,2),X_departure_orbit(1,3),'o','MarkerEdgeColor','none','MarkerFaceColor','green','MarkerSize',8);
StartRad_ArrOrbit = plot3(X_arrival_orbit(1,1),X_arrival_orbit(1,2),X_arrival_orbit(1,3),'o','MarkerEdgeColor','none','MarkerFaceColor','red','MarkerSize',8);
transf_orbit = plot3(X_transf(:,1),X_transf(:,2),X_transf(:,3),'y');
departure_orbit = plot3(X_departure_orbit(:,1),X_departure_orbit(:,2),X_departure_orbit(:,3),'g');
arrival_orbit = plot3(X_arrival_orbit(:,1),X_arrival_orbit(:,2),X_arrival_orbit(:,3),'r');
RadMan_1 = plot3(R_maneouvre_1(1,f),R_maneouvre_1(2,f),R_maneouvre_1(3,f),'o','MarkerEdgeColor','none','MarkerFaceColor','blue','MarkerSize',8);
RadMan_2 = plot3(R_maneouvre_2(1,f),R_maneouvre_2(2,f),R_maneouvre_2(3,f),'o','MarkerEdgeColor','none','MarkerFaceColor','y','MarkerSize',8);
Man2Rad_DepOrbit = plot3(X_dep_t02Tman2(end,1),X_dep_t02Tman2(end,2),X_dep_t02Tman2(end,3),'o','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',8);
Man1Rad_ArrOrbit = plot3(X_arr_t02Tman1(end,1),X_arr_t02Tman1(end,2),X_arr_t02Tman1(end,3),'o','MarkerEdgeColor','none','MarkerFaceColor',[0.4,0.1,0.3],'MarkerSize',8);
xlabel('x')
ylabel('y')
zlabel('z')
title('Orbit Representation','FontSize',13)
legend([transf_orbit,departure_orbit,arrival_orbit,StartRad_DepOrbit,StartRad_ArrOrbit,RadMan_1,RadMan_2,Man2Rad_DepOrbit,Man1Rad_ArrOrbit],{'transfer orbit','departure orbit','arrival orbit',...
    'start position for the departure orbit','start position for the arrival orbit','first maneouvre point','second maneouvre point','position of pl_1 @second maneouvre time','position of pl_2 @first maneouvre time'})

I = imread('Sun.png'); 
Sun_radius = astroConstants(3);
[x_sun, y_sun, z_sun] = ellipsoid (0,0,0,Sun_radius,Sun_radius,Sun_radius);
Sun = surf(20*x_sun, 20*y_sun, 20*z_sun,'Edgecolor', 'none');
set(Sun,'FaceColor','texturemap','Cdata',I)
% set(gca,'Color','none','visible','off')
axis equal

hold off


elapsed_time = toc;  
fprintf('ELAPSED TIME: %g [s]\n\n\n',elapsed_time)