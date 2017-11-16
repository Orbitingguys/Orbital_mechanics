r_geo = [10000 15000 13000];
v_geo = [1 1 1];
mi = 398600;

[orbital_parameters_1] = geo2kep(r_geo,v_geo,mi);
[r_geo_f,v_geo_f] = kep2geo (orbital_parameters_1,mi);