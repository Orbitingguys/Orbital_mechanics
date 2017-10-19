%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rot_matrix.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function che calcola la matrice di rotazione il sistema perifocale e
% quello geocentrico.
% 
% Input:  - incli    [rad]            inclinazione desiderata;
% 
%         - omegone: [rad]            RAAN voluta:
% 
%         - omeghino [rad]            argomento del pericentro da ottenere;
% 
% Output: - R        [adimensionale]  matrice di rotazione definita dai 3
%
%                                     angoli in ingresso;
% 
% Autori:    Adriano Filippo Inno | Patrizia Impieri | Piergiorgio Frulla
% Matricole:      (828084)              (825837)           (827798)
% Prova finale di Analisi di Missioni spaziali, gruppo A19 
%-------------------------------------------------------------------------%

function   [R] = rot_matrix_geo2pf(INCLI, RAAN, PA)

% if abs(omegone) < 2*pi && abs(omeghino) < 2*pi  && abs(incli) < 2*pi 
% rot omegone (attorno a z)
R_omegone = [  cos(RAAN)    sin(RAAN)     0  ;
              -sin(RAAN)    cos(RAAN)     0  ;
                     0               0          1 ];
                 
% rot incli (attorno a x')                
R_incli =     [      1       0             0       ;
                     0   cos(INCLI)    sin(INCLI)  ;
                     0  -sin(INCLI)    cos(INCLI) ];

% rot omeghino (attorno a z'')
R_omeghino = [  cos(PA)   sin(PA)   0  ;
               -sin(PA)   cos(PA)   0  ;
                    0                0          1] ;
                
% else
%     R_omegone = [  cosd(omegone)    sind(omegone)     0  ;
%               -sind(omegone)    cosd(omegone)     0  ;
%                      0               0          1 ];
%                  
% % rot incli (attorno a x')                
% R_incli =     [      1       0             0       ;
%                      0   cosd(incli)    sind(incli)  ;
%                      0  -sind(incli)    cosd(incli) ];
% 
% % rot omeghino (attorno a z'')
% R_omeghino = [  cosd(omeghino)   sind(omeghino)   0  ;
%                -sind(omeghino)   cosd(omeghino)   0  ;
%                     0                0          1] ;

R = R_omeghino*R_incli*R_omegone;


end