function  [r_geo,v_geo] = kep2geo (orbital_parameters,mi)

a = orbital_parameters.a;
ecc = orbital_parameters.ecc;
RAAN = orbital_parameters.RAAN;
PA = orbital_parameters.PA;
INCLI = orbital_parameters.INCLI;
theta = orbital_parameters.theta;





p = a*(1-ecc^2);
% norm_h = sqrt(p*mi);
norm_r = p./(1+ecc*cos(theta));
r_pf = norm_r.*[cos(theta); sin(theta); 0];
[R] = geo2perif_rot(RAAN,PA,INCLI);
r_geo = R' * r_pf;
v_solidal = sqrt(mi/p).*[ecc*sin(theta) 1+ecc*cos(theta) 0];

v_pf = [v_solidal(1)*cos(theta)-v_solidal(2)*sin(theta); v_solidal(1)*sin(theta)+v_solidal(2)*cos(theta); 0];
v_geo = R' * v_pf;

% h_pf = [0 0 norm_h];
% h_geo = R' * h_pf;
% e_pf = [ecc 0 0];
% e_geo = R' * e_pf;
