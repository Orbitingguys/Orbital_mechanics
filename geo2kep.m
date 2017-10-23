function [orbital_parameters] = geo2kep(r_geo,v_geo,mi)

I_geo = [1, 0, 0];            % [km] versore asse X_geo
J_geo = [0, 1, 0];            % [km] versore asse Y_geo
K_geo = [0, 0, 1];
norm_r = norm(r_geo);
norm_v = norm(v_geo);
a = 1/(2/norm_r-norm_v.^2/mi);
h = cross(r_geo,v_geo);
norm_h = norm(h);
e = (cross(v_geo,h) - mi.*r_geo./norm_r)./mi;
ecc = norm(e);
% if ecc <= 0.0001
%     fprintf('----Geometria orbita---- \n\n')
%     fprintf('orbita circolare \n\n\n');
% elseif ecc > 0.0001 && ecc <= 0.9999
%     fprintf('----Geometria orbita----: \n\n')
%     fprintf('orbita ellittica \n\n\n');
% elseif ecc > 0.9999 && ecc <= 1.0001
%     fprintf('----Geometria orbita----: \n\n')
%     fprintf('orbita parabolica \n\n\n');
% else
%     fprintf('----Geometria orbita----: \n\n')
%     fprintf('orbita iperbolica \n\n\n');
% end
INCLI = acos(dot(h,K_geo)/norm_h);
N = cross(K_geo,h)/norm(cross(K_geo,h));
if dot(N,J_geo)>=0 
    RAAN = acos(dot(I_geo,N));
else
    RAAN = (2*pi-acos(dot(I_geo,N)));
end

if dot(e,K_geo)>=0 
    PA = acos(dot(N,e)/ecc);
else
    PA = (2*pi-acos(dot(N,e)/ecc));
end

if dot(v_geo,r_geo)>=0 
    theta = acos(dot(r_geo,e)/(norm_r*ecc));
else
    theta = (2*pi-acos(dot(r_geo,e)/(norm_r*ecc)));
end
orbital_parameters.a = a;
orbital_parameters.ecc = ecc;
orbital_parameters.INCLI = INCLI;
orbital_parameters.RAAN = RAAN;
orbital_parameters.PA = PA;
orbital_parameters.theta = theta;