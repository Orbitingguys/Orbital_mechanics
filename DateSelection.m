function [DateOption] = DateSelection(DateSetup)
%% DateSetup
%                       .default:       -   If .default is 1 (True) the dates are directly taken
%                                           from the data given for each body.
%
%                                       -   If .default is 0 (False), DEPENDING ON
%                                        DefaultOwn_Date.TypeMission, date must be defined by the
%                                        varables:
%                                                   .date_d_min:    First possible date for departure
%                                                   .date_d_max:    Last possible date for departure
%                                                   .date_a_min:    First possible date for arrival
%                                                   .date_a_max:    Lats possible date for arrival
%
%
%                                           FOR THE CASE DateSetup.default==1 NO MORE VARIABLES
%                                           ARE NEEDED (EXCEPT DefaultOwn_Date.jDatebody).


%%

switch DateSetup.default
    
    case 1
        
        [~,~,DateOption]=body_time_data(DateSetup.jbody);
        
    case 0
        
        DateOption.TypeMission=DateSetup.TypeMission;
        DateOption.date_d_min = DateSetup.date_d_min;
        
        switch DateOption.TypeMission
    
    
            case 1
        
                DateOption.date_d_max=DateSetup.date_d_max;
                DateOption.date_a_min=DateSetup.date_a_min;
                DateOption.date_a_max=DateSetup.date_a_max;
        
            case 2
        
                DateOption.date_d_max=DateSetup.date_d_max;
        
            case 3
        
                DateOption.date_a_min=DateSetup.date_a_min;
                DateOption.date_a_max=DateSetup.date_a_max;
        
            case 4
       
                %Nothing more is required.
        end
       
end