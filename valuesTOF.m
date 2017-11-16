function [vt] = valuesTOF(T_d,T_a)
%% TIME OF FLIGHT LINES

% function that computes the lines that defines the TOF for the contour plot
% the order of the minimum and maximun TOF are computed in order to divide 
% properly the interval of values

%% MINIMUM OF TOF

minTOF = T_a(1)-T_d(end); % this value can be less or equal to 0!

if minTOF < 0 % the first contour line will have a negative value
    
    Order_minTOF = floor(log10(-minTOF)); % computing the order of the minimum TOF
    
elseif minTOF == 0 % the first contour line will start at 0
    
    Order_minTOF = 1;
    
else % the first contour line will have a positive value 
    
    Order_minTOF=floor(log10(minTOF));
    
end

%% MAXIMUM OF TOF

maxTOF = T_a(end)-T_d(1); % this value is always be positive

Order_maxTOF = floor(log10(maxTOF)); % computing the order of the maximum TOF

%% COMPUTING TOF LINES

startTOF = floor(minTOF/10^(Order_minTOF-1))*10^(Order_minTOF-1); % value of the first contour line that will be contour plotted
endTOF = ceil(maxTOF/10^(Order_maxTOF-1))*10^(Order_maxTOF-1);      % value of the last contour line that will be contour plotted
n = 10;                                                             % number of lines we want to apear in the contour plot

vt = round(linspace(startTOF,endTOF,n));                            % vector needed for the contour rounding the days to the unity.


