function [orbital_parameters, delta_v] = half_hyp(v_i_inc_or_out, v_f_out_or_inc, R_body, V_body, mi_body, orbital_parameters_park, orbital_parameters_body, option)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% half_hyp.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
%                                                                                                                    %     
% INPUTS:    v_i_inc_or_out [3]         velocity of entering the hyperbola in the heliocentric reference frame [Km/s]%     
%            v_f_inc_or_out [3]         velocity of exiting the hyperbola in the heliocentric reference frame [Km/s] % 
%            R_body [3]                 position of the planet in the heliocentric reference frame [Km]              % 
%            mi_body [1]                gravitational constant of the body [Km^3/s^2]                                %   
%            option                     - = 1: from a parking orbit to an interplanetary orbit                       %    
%                                       - = 2: from an interplanetary orbit to a parking orbit (rendez-vouz)         %                                                  
%                                                                                                                    %
% OUTPUTS:   orbital_parameters{6}      - .a, major semiaxis [km]                                                    %    
%                                       - .ecc, eccentricity [/]                                                     % 
%                                       - .RAAN, right ascension of ascending node [rad]                             %    
%                                       - .PA, pericenter argument [rad]                                             %    
%                                       - .INCLI, inclination of the orbit [rad]                                     %    
%                                       - .theta_inf, the infinity theta of the hyperbola [rad]                      %    
%                                       - .incoming, if equal to 1 the cycle computes the incoming                   %                                                 
%                                         hyperbola (infinity velocity minus) and viceversa if it is equal to 0      %                                                                                                                    %   
%                                                                                                                    %
% AUTORS:   ADRIANO FILIPPO INNO               PIERGIORGIO FRULLA          PALLAR?S CHAMORRO FRANCISCO DE ASIS       %                                                              %   
%                                                                                                                    %
%           PROJECT-ORBITAL MECHANICS                  2017/2018                   GRUOP NUMBER 19                   %                                 
%                                                                                                                    %       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = 6.67259e-20;
msun = 1.988919445342813e+030;
mi_sun = msun*G;

% r_vers = -R_body/norm(R_body);   % directed towards the Sun
% t_vers = V_body/norm(V_body);
% k_vers = cross(-r_vers,t_vers)/norm(cross(-r_vers,t_vers));
% theta_vers = cross(-r_vers,k_vers);

a_body = orbital_parameters_body.a;
ecc_body = orbital_parameters_body.ecc;
RAAN_body = orbital_parameters_body.RAAN;
PA_body = orbital_parameters_body.PA;
INCLI_body = orbital_parameters_body.INCLI;
theta_body = orbital_parameters_body.theta;
p_body = a_body*(1-ecc_body^2);

a_park = orbital_parameters_park.a;
ecc_park = orbital_parameters_park.ecc;
RAAN_park = orbital_parameters_park.RAAN;
PA_park = orbital_parameters_park.PA;
INCLI_park = orbital_parameters_park.INCLI;
theta_start_park = orbital_parameters_park.theta;
p_park = a_park*(1-ecc_park^2);

if option == 1
    v_i_out = v_i_inc_or_out;
    r_p = p_park/(1+ecc_park*cos(theta_start_park));
    v_f_out = v_f_out_or_inc;
    v_inf_plus = v_f_out - V_body;
    a = -mi_body/(norm(v_inf_plus))^2;
    ecc = 1 + (r_p*(norm(v_inf_plus))^2)/mi_body;
    theta_inf = atan2(-1/ecc,sqrt(ecc^2-1)/ecc);    % from -pi to pi
    orbital_parameters.incoming = false;
    
    
elseif option == 2
    v_i_inc = v_i_inc_or_out;
    v_f_inc = v_f_out_or_inc;
    r_p = p_park/(1+ecc_park*cos(theta_start_park));
    v_inf_minus = v_i_inc - V_body;
    a = -mi_body/(norm(v_inf_minus))^2;
    ecc = 1 + (r_p*(norm(v_inf_minus))^2)/mi_body;
    theta_inf = atan2(-1/ecc,sqrt(ecc^2-1)/ecc);
    orbital_parameters.incoming = true;
end

orbital_parameters.a = a;
orbital_parameters.ecc = ecc;
orbital_parameters.RAAN = 0;
orbital_parameters.PA = 0;
orbital_parameters.INCLI = 0;
orbital_parameters.theta_inf = theta_inf;

delta_v = 0;