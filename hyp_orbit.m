function [r_body2sc,r_pf] = hyp_orbit(orbital_parameters,mi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyp_orbit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
%                                                                                                                %     
% INPUTS:   orbital_parameters[7]       - .a, major semiaxis [km]                                                %     
%                                       - .ecc, eccentricity [/]                                                 % 
%                                       - .RAAN, right ascension of ascending node [rad]                         %    
%                                       - .PA, pericenter argument [rad]                                         %    
%                                       - .INCLI, inclination of the orbit [rad]                                 %    
%                                       - .theta_inf, the infinity theta of the hyperbola [rad]                  %    
%                                       - .incoming, if equal to 1 the cycle computes the incoming               %                                                 
%                                         hyperbola (infinity velocity minus) and viceversa if it is equal to 0  %    
%                                                                                                                %
% OUTPUTS:  r_body2sc[3]                  radius in the inertial frame of the body that is being considered [km] %    
%           r_pf[3]                       radius in the perifocal frame                                          %   
%                                                                                                                %
% AUTORS:   ADRIANO FILIPPO INNO               PIERGIORGIO FRULLA          PALLAR?S CHAMORRO FRANCISCO DE ASIS   %                                                              %   
%                                                                                                                %
%           PROJECT-ORBITAL MECHANICS                  2017/2018                   GRUOP NUMBER 19               %                                 
%                                                                                                                %       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                    

a = orbital_parameters.a;
ecc = orbital_parameters.ecc;
RAAN = orbital_parameters.RAAN;
PA = orbital_parameters.PA;
INCLI = orbital_parameters.INCLI;
theta_inf = orbital_parameters.theta_inf;

p = a*(1-ecc^2);
[R] = body2perif_rot(INCLI, RAAN, PA);

switch nargin
    case 3
%         phi = atan2(R_sun2body(1),R_sun2body(2));
    case 2 
%         phi = 0;
end


if a > 0
    
    error('the semimajor axes must be negative, do you mean a_bar?')
    
end

if ecc <= 0
    
    error('the eccentricity of a hyperbola must be at least 1')
    
end

% if theta_inf >= 2*pi 
%     
%     error('theta infinity maybe is in degrees')
%     
% elseif theta_inf <= 0 
%     
%     error('theta infinity is less that zero')
%     
% end


if orbital_parameters.incoming == 1
    
    theta = linspace(theta_inf,0,10000);
    
else
    
    theta = linspace(0,theta_inf,10000);
    
end

% pre_allocation of the memory
r_pf = zeros(3,length(theta));
r_body2sc = zeros(3,length(theta));
v_solidal = zeros(3,length(theta));
v_pf = zeros(3,length(theta));
v_body2sc = zeros(3,length(theta));

for i = 1:length(theta)
    
    norm_r = p/(1+ecc*cos(theta(i))); 
    r_pf(:,i) = norm_r.*[cos(theta(i)); sin(theta(i)); 0];                    % radius in the perifocal frame
    r_body2sc(:,i) = R'*r_pf(:,i);                                            % radius in the inertial frame of the body
    v_solidal(:,i) = sqrt(mi/p).*[ecc*sin(theta(i)) 1+ecc*cos(theta(i)) 0];   % velocity vector solidal to the s/c motion
    v_pf(:,i) = [v_solidal(1,i)*cos(theta(i))-v_solidal(2,i)*sin(theta(i));   % velocity vector in the perifocal frame
                 v_solidal(1,i)*sin(theta(i))+v_solidal(2,i)*cos(theta(i)); 
                 0                                                       ];
    v_body2sc(:,i) = R'*v_pf(:,i);                                            % velocity vector in the inertial frame of the body
%     r_helio(:,i) = norm(r_body2sc(:,i)).*[cos(phi); sin(phi); 0]; 
%     v_helio(:,i) = norm(v_body2sc(:,i)).*[cos(phi); sin(phi); 0];
    
end





    
    