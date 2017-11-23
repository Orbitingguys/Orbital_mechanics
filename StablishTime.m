function [TimeOption] = StablishTime(DateOption)



% DateOption     .TypeMission:            FOR EVERY CASE, THE FIRST POSSIBLE DATE OF DEPARTURE .date_d_min
%                                         MUST BE ADDED FOR COMPUTING THE INCREMENTS.
%
%                                          -   if .TypeMission is 1 it means that you have a window of
%                                              launch so you need to define .date_d_max and also you have a
%                                              request for the arrival date which add a constrain that has
%                                              to be saved in matlab as .date_a_min and .date_a_max.
%
%                                          -   if .TypeMission is 2 you still have the launch window
%                                              but you don't deserve a constrain for the arrival date, so
%                                              .date_d_max will be needed.
%
%                                          -   if .TypeMission is 3 you have to respect only an arrival
%                                              date, so .date_a_min and .date_a_max will needed.
%
%                                          -   if .TypeMission is 4 you haven't any time constrain.
%   
%                .date_d_max:              -   Must be added for .TypeMission equal to 1 or 2
%                                              for the other cases will be computed as the period of the
%                                              first orbit.
%                                       
%                .date_a_min, date_a_max:  -   Must be added for .TypeMission equal to 1 or 3
%                                              for the other cases will be computed as 3 times of the
%                                              period of the second orbit.
%
%                WHEN A DATE IS NOT NEEDED JUST DON'T ADD IT AS A DATA!

TimeOption.TypeMission = DateOption.TypeMission;
date_d_min = DateOption.date_d_min; 

date0 = date_d_min;                         % the date of reference neded to compute the increments of time.

TimeOption.mjd2000_0 = date2mjd2000(date0); % mjd2000 of date0
mjd2000_0 = TimeOption.mjd2000_0;

switch TimeOption.TypeMission
    
    case 1
        
        date_d_max = DateOption.date_d_max;
        date_a_min = DateOption.date_a_min;
        date_a_max = DateOption.date_a_max;
        
        % calculating the mjd2000 of the dates.
        mjd2000_d_max = date2mjd2000(date_d_max);
        mjd2000_a_min = date2mjd2000(date_a_min);
        mjd2000_a_max = date2mjd2000(date_a_max);
        
        % finally, we compute the time as the increments of days (mjd2000) since date0.
        TimeOption.t_i_max = (mjd2000_d_max-mjd2000_0)*3600*24;
        TimeOption.t_f_min = (mjd2000_a_min-mjd2000_0)*3600*24;
        TimeOption.t_f_max = (mjd2000_a_max-mjd2000_0)*3600*24;
        
    case 2
        
        date_d_max = DateOption.date_d_max;
        
        % calculating the mjd2000 of the dates.
        mjd2000_d_max = date2mjd2000(date_d_max);
        
        % finally, we compute the time as the increments of days (mjd2000) since date0.
        TimeOption.t_i_max = (mjd2000_d_max-mjd2000_0)*3600*24;
        
    case 3
        
        date_a_min = DateOption.date_a_min;
        date_a_max = DateOption.date_a_max;
        
        % calculating the mjd2000 of the dates.
        mjd2000_a_min = date2mjd2000(date_a_min);
        mjd2000_a_max = date2mjd2000(date_a_max);
        
        % finally, we compute the time as the increments of days (mjd2000) since date0.
        TimeOption.t_f_min = (mjd2000_a_min-mjd2000_0)*3600*24;
        TimeOption.t_f_max = (mjd2000_a_max-mjd2000_0)*3600*24;
        
    case 4
       
        %Nothing more is required.
       
end
    
 