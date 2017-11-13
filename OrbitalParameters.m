function [orbital_parameters]=OrbitalParameters(mjd2000,body)
if body<9
    
    [kep,~] = uplanet(mjd2000,body);
    orbital_parameters.a=kep(1);
    orbital_parameters.ecc=kep(2);
    orbital_parameters.INCLI=kep(3);
    orbital_parameters.RAAN=kep(4);
    orbital_parameters.PA=kep(5);
    orbital_parameters.theta=kep(6);
    
else
    
    [kep,~,~] = ephNEO(mjd2000,body);
    orbital_parameters.a=kep(1);
    orbital_parameters.ecc=kep(2);
    orbital_parameters.INCLI=kep(3);
    orbital_parameters.RAAN=kep(4);
    orbital_parameters.PA=kep(5);
    orbital_parameters.theta=kep(6);
    
end