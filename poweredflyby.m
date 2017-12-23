function [orbital_parameters1, orbital_parameters2, delta_v]= poweredflyby (V_i, V_f,orbital_parameters_planet,mi_sun,mi_planet)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hyp_orbit.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
%                                                                                                                %     
% INPUTS:   orbital_parameters_planet[7]    - .a, major semiaxis of planet @sun [km]                             %     
%                                           - .ecc, eccentricity of planet @sun [/]                              % 
%                                           - .RAAN, right ascension of ascending node of planet @sun [rad]      %    
%                                           - .PA, pericenter argument of planet @sun [rad]                      %    
%                                           - .INCLI, inclination of the orbit of planet @sun [rad]              %    
%                                           - .theta, true anomally of planet @sun [rad]                         %    
%            V_i:       Incoming velocity in the heliocentric ecliptic frame                                     %
%            V_f:       Outgoing velocity in the heliocentric ecliptic frame                                     %
%            mi_sun:    Gravitational parameter of the sun                                                       %
%            mi_planet: Gravitational parameter of the planet                                                    %
%                                                                                                                %
% OUTPUTS:  orbital_parameters1: Keplerean parameters of the initial hyperbola related to the heliocentric       %
%                                ecliptic frame                                                                  %
%           orbital_parameters2: Keplerean parameters of the initial hyperbola related to the heliocentric       %
%                                ecliptic frame                                                                  %
%           delta_v:             delta v powered by the thrusters at the pericenter of the hyperbolas            %
%                                                                                                                %
% AUTORS:   ADRIANO FILIPPO INNO               PIERGIORGIO FRULLA          PALLARÉS CHAMORRO FRANCISCO DE ASIS   %   
%                                                                                                                %
%           PROJECT-ORBITAL MECHANICS                  2017/2018                   GRUOP NUMBER 19               %                                 
%                                                                                                                %       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                    

[Rp,V_p] = kep2geo (orbital_parameters_planet,mi_sun);
vinf_i=V_i-V_p;
vinf_1=norm(vinf_i);
gammai=atan2(vinf_i(2),vinf_i(1));
vinf_f=V_f-V_p;
vinf_2=norm(vinf_f);
gammaf=atan2(vinf_f(2),vinf_f(1));
delta=abs(gammaf-gammai);
if delta>pi
    delta=2*pi-delta;
end
x_0=[1.1 1.1 pi/3 pi/3 1000]';
X=fsolve(@(x)hyp_system(x,vinf_2,vinf_1,delta,mi_planet),x_0);
ecc1=X(1);
ecc2=X(2);
delta1=X(3);
delta2=X(4);
rp=X(5);
n=cross(vinf_i,vinf_f);n=n/norm(n);
ln=cross([0;0;1],n);ln=ln/norm(ln);
%Complete it with the parallel vector case
jn=cross(n,ln);
INCLI=acos(n(3));
if INCLI==0 || INCLI==pi
    RAAN=0;
    ln=[1; 0; 0];
else
    if ln(2)>=0
        RAAN=acos(ln(1));
    else
        RAAN=2*pi-acos(ln(1));
    end
end
%We define alfa as the angle between the line of nodes and the correspondant velocity at infinity on
%a hyperbolic trajectory
if vinf_i'*jn>=0
    alfa1=acos(dot(vinf_i,ln)/vinf_1);
else
    alfa1=-acos(dot(vinf_i,ln)/vinf_1);
end
if vinf_f'*jn>=0
    alfa2=acos(dot(vinf_f,ln)/vinf_2);
else
    alfa2=-acos(dot(vinf_f,ln)/vinf_2);
end
%We calculate the true anomallies at infinity points
thetainf_1=-acos(-1/ecc1);
thetainf_2=acos(-1/ecc2);
PA1=alfa1-thetainf_1-pi;
PA2=alfa2-thetainf_2;
rad2deg(PA1-PA2)
a1=rp/(1-ecc1);
a2=rp/(1-ecc2);
orbital_parameters1.a=a1;
orbital_parameters1.ecc=ecc1;
orbital_parameters1.RAAN=RAAN;
orbital_parameters1.INCLI=INCLI;
orbital_parameters1.PA=PA1;
orbital_parameters2.a=a2;
orbital_parameters2.ecc=ecc2;
orbital_parameters2.RAAN=RAAN;
orbital_parameters2.INCLI=INCLI;
orbital_parameters2.PA=PA1;
vp_1=sqrt(vinf_1^2+2*mi_planet/rp);
vp_2=sqrt(vinf_2^2+2*mi_planet/rp);
delta_v=vp_2-vp_1;
