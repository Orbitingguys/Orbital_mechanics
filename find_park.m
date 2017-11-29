function [orbital_parameters] = find_park(R_sun2body,V_body,v_inc,mi_body,R_body,SOI)

v_inf_minus = v_inc - V_body;
v_inf_minus_norm = norm(v_inf_minus);

h_vers = cross(R_sun2body,v_inf_minus)/norm(cross(R_sun2body,v_inf_minus));

if h_vers(3) > 0
    INCLI = acos(dot(h_vers,[0 0 1]));
else 
    INCLI = pi -acos(dot(h_vers,[0 0 1]));
end
N = cross([0 0 1],h_vers)/norm(cross([0 0 1],h_vers));
if dot(N,[0 1 0]) >= 0 
    RAAN = acos(dot([1 0 0],N));
else
    RAAN = (2*pi-acos(dot([1 0 0],N)));
end


delta_v = @(rp) sqrt((2*mi_body+1+rp.*v_inf_minus_norm^2)./(2.*rp))-sqrt(mi_body./rp);
r_p = ga(delta_v,1,[],[],[],[],R_body,SOI);

rp = R_body:0.5:SOI;
plot(rp,delta_v(rp))

orbital_parameters.ecc = 0;
orbital_parameters.a = r_p;
orbital_parameters.INCLI = INCLI;
orbital_parameters.RAAN = RAAN;

