function [printed_value,ERROR]=PORKCHOP_PROCEDURE(DateSetup,ibody,jbody)
%% INPUTS
%                   -   ibody:      Body from the Solar System from which the mission departs
%
%                   -   jbody:      Body from the Solar System to which the mission arrives
%
%                   -   DateSetup:  Starting defining DateSetup.default
%
%                       .default:       -   If .default is 1 (True) the dates are directly taken
%                                           from the data given for each body.
%
%                                       -   If .default is 0 (False), DEPENDING ON
%                                        DefaultOwn_Date.TypeMission, date must be defined by the
%                                        varables:
%                                                   .date_d_min:    First possible date for departure
%                                                   .date_d_max:    Last possible date for departure
%                                                   .date_a_min:    First possible date for arrival
%                                                   .date_a_max:    Last possible date for arrival
%
%
%                                           FOR THE CASE DateSetup.default==1 NO MORE VARIABLES
%                                           ARE NEEDED.
%%

G=6.67259e-20;
msun=1.988919445342813e+030;
mi=msun*G;

DateSetup.ibody=ibody;
DateSetup.jbody=jbody;
[DateOption] = DateSelection(DateSetup);
T0=date2mjd2000(DateOption.date_d_min)*3600*24;
[TimeOption]=StablishTime(DateOption);
[orbital_parameters_1]=OrbitalParameters(TimeOption.mjd2000_0,ibody);
[orbital_parameters_2]=OrbitalParameters(TimeOption.mjd2000_0,jbody);
[printed_value,ERROR] = pork_chopFORWHILEedit2(orbital_parameters_1,orbital_parameters_2,TimeOption,mi,T0);