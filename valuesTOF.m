function [vt]=valuesTOF(T_d,T_a)

%TIME OF FLIGHT LINES
%We wanto to create the contour plot with the lines that defines the TOF of
%the points, for that we need to know the order of the minimum and maximun
%TOF in order to divide properly de interval ovalues
minTOF=T_a(1)-T_d(end);
%This value can be less or equal to 0, a fact that can bring us problems
%after
if minTOF<0 %Then, the first contour line will have a negative value
    Order_minTOF=floor(log10(-minTOF));
elseif minTOF==0 %Then, the first contour line will start at 0
    Order_minTOF=1;
else %Then, the first contour line will have a positive value
    Order_minTOF=floor(log10(minTOF));
end
maxTOF=T_a(end)-T_d(1);
%In the other hand, this value must always be positive, otherwise the
%analysis have no sense.
Order_maxTOF=floor(log10(maxTOF));
%We will now define the values of the first and last contour lines we want
%the contour plot to represent, for that we floor the minimum TOF / ceil
%the maximun TOF to the lower / upper ten
startTOF=floor(minTOF/10^(Order_minTOF-1))*10^(Order_minTOF-1);
endTOF=ceil(maxTOF/10^(Order_maxTOF-1))*10^(Order_maxTOF-1);
%Then, we decide how many lines we want to apear in the contour plot, for
%example 10
n=10;
%and finally, we create the vector needed as the contour doc says, rounding
%the days to the unity.
vt=round(linspace(startTOF,endTOF,n));