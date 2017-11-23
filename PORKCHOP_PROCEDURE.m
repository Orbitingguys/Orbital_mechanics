function [printed_value,ERROR] = PORKCHOP_PROCEDURE(DateSetup)
%% INPUTS
%                   -   ibody:      Body from the Solar System from which the mission departs
%
%                   -   jbody:      Body from the Solar System to which the mission arrives
%
%                   -   DateSetup:  Starting defining DateSetup.default
%
%                       .default:       -   If .default is 1 (True) the dates are directly taken
%                                           from the data given for each body.
%
%                                       -   If .default is 0 (False), DEPENDING ON
%                                           DefaultOwn_Date.TypeMission, date must be defined by the
%                                           varables:
%                                                   .date_d_min:    First possible date for departure
%                                                   .date_d_max:    Last possible date for departure
%                                                   .date_a_min:    First possible date for arrival
%                                                   .date_a_max:    Last possible date for arrival
%
%
%                                       FOR THE CASE DateSetup.default==1 NO MORE VARIABLES
%                                       ARE NEEDED.
%
%                   REMEMBER TO ADD TO MATLAB PATH THE FOLDER 'TimeConversion'  

%%
G = 6.67259e-20;
msun = 1.988919445342813e+030;
mi = msun*G;
% 
% ibody = DateSetup.ibody;
% jbody = DateSetup.jbody;
[DateOption] = DateSelection(DateSetup);
T0 = date2mjd2000(DateOption.date_d_min)*3600*24;
[TimeOption] = StablishTime(DateOption);
[orbital_parameters_1] = OrbitalParameters(TimeOption.mjd2000_0,DateSetup.ibody);
[orbital_parameters_2] = OrbitalParameters(TimeOption.mjd2000_0,DateSetup.jbody);
m = round((TimeOption.t_i_max)/86400);
n = round((TimeOption.t_f_max-TimeOption.t_f_min)/86400);

fprintf('Mission from %s',char(DateOption.target_i))
fprintf(' to %s',char(DateOption.target_j))
fprintf(' is being analyzed \n\n')

if m+n < 600
    
    [printed_value,ERROR] = pork_chopFORWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi,T0,m,n);
    
elseif m+n >= 15000
    
    while m+n>15000
        m = round(m/10);
        n = round(n/10);
    end
    
    if m < 5
       m = m + 10;
    end

    [printed_value,ERROR] = pork_chopPARPOLFORWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi,T0,m,n);
    
elseif m+n > 2500 && m+n < 15000
    
    while m+n > 2500
        m = round(m/2);
        n = round(n/2);
    end

    [printed_value,ERROR] = pork_chopPARPOLFORWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi,T0,m,n);
    
else

    [printed_value,ERROR] = pork_chopPARPOLFORWHILE(orbital_parameters_1,orbital_parameters_2,TimeOption,mi,T0,m,n);

end


