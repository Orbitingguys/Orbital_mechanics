function HeaderExample
%function [Xd] = Gaussplanetary (t,X,T0)

% Gaussplanetary.m - This function returns the derivative of the vector of six osculating 
%           parameters, at a point defined by those parameters on a perturbed problem.
%
% PROTOTYPE:
%   [Xd] = Gaussplanetary (t,X,T0)
%
% INPUT:
%   t[1]    Time since T0. [s].
%   X[6]    Vector containing the six osculating parameters, by order: Specific
%           Angular Momentum (h) [km^2/s], Eccentricity (ecc) [-], True Anomally (theta) [deg],Right
%           Ascension of the Ascending Node (RA) [deg], Inclination (i) [deg], Argument of Periapse 
%           (om) [deg].
%   T0[1]   Time since 01-01-2000, 12:00 noon. [s].
%   
% OUTPUT:
%   Xd[6]   Vector contining the derivative of the six osculating parameters, by order: Specific
%           Angular Momentum (h), Eccentricity (ecc), True Anomally (theta),Right Ascension of the 
%           Ascending Node (RA), Inclination (i), Argument of Periapse (om).
%
% CALLED FUNCTIONS:
%   astroConstants.m
%   J2pert.m
%   moonpert.m
%
% REFERECES:
%
% AUTHORS:
%   Francisco Pallarés Chamorro, 29/01/2018, MATLAB, Gausplanetary.m
%---------------------------------------------------------------------------------------------------