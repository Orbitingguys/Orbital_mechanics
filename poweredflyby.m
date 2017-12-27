function [orbital_parameters1, orbital_parameters2, delta_v]= poweredflyby (V_i, V_f,orbital_parameters_planet,mi_sun,mi_planet,RPmin)


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
vinf_f=V_f-V_p;
vinf_2=norm(vinf_f);
delta=atan2(norm(cross(vinf_i,vinf_f)),dot(vinf_i,vinf_f));
if delta>pi
    delta=2*pi-delta;
end
options = optimset('Display','off','TolFun',1e-16);
rp=fsolve(@(rp)hyp_system_simple(rp,vinf_2,vinf_1,delta,mi_planet),RPmin,options);
if rp>RPmin
    %We can calculate directly the velocities at the periapse
    vp_1=sqrt(vinf_1^2+2*mi_planet/rp);
    vp_2=sqrt(vinf_2^2+2*mi_planet/rp);
    %And the delta_v needed as well
    delta_v=vp_2-vp_1;
    %We proceed now to compute the parameters of the orbits
    ecc1=1+rp*vinf_1^2/mi_planet;
    ecc2=1+rp*vinf_2^2/mi_planet;
    h1=rp*vp_1;
    h2=rp*vp_2; 
    delta1=2*asin(1/ecc1);
    delta2=2*asin(1/ecc2);
    n=cross(vinf_i,vinf_f);n=n/norm(n);
    h1_v=h1.*n;
    h2_v=h2.*n;
    ecc1_v=cross(vinf_i,h1_v)/mi_planet+vinf_i/vinf_1;
    ecc1_vnorm=ecc1_v./norm(ecc1_v);
    ecc2_v=cross(vinf_f,h2_v)/mi_planet-vinf_f/vinf_2;
    ecc2_vnorm=ecc2_v./norm(ecc2_v);
    %Computation of the INCLI using the normal vectors to the ecliptic plane an the trajectory plane
    ln=cross([0;0;1],n);ln=ln/norm(ln);
    INCLI=acos(n(3));
    %Computation of the RAAN using the vectors of the line of nodes and the gamma point
    if INCLI==0 || INCLI==pi %In case of trajectory laying on the ecliptic plane
        RAAN=0;
        ln=[1; 0; 0];
    else
        if ln(2)>=0
            RAAN=acos(ln(1));
        else
            RAAN=2*pi-acos(ln(1));
        end
    end
    %Computation of the PA using directly the eccentricity vectors and the line of nodes
    if ecc1_vnorm(3)>=0
        PA1=acos(dot(ln,ecc1_vnorm));
    else
        PA1=2*pi-acos(dot(ln,ecc1_vnorm));
    end
    if ecc2_vnorm(3)>=0
        PA2=acos(dot(ln,ecc2_vnorm));
    else
        PA2=2*pi-acos(dot(ln,ecc2_vnorm));
    end
    %Analysis of the difference between both computed PA (They should be equal)
    ERROR=rad2deg(PA1-PA2); %[º]
    if ERROR>=0.1
        rp
        pa1=rad2deg(PA1)
        pa2=rad2deg(PA2)
        ERROR=ERROR
    end
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
else
    orbital_parameters1=NaN;
    orbital_parameters2=NaN;
    delta_v=NaN;
end