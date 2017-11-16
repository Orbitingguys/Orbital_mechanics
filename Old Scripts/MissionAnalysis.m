function [Dvtotmin,T_min_d,T_min_a,T_d,T_a,Tof]=MissionAnalysis(date_min_d,date_max_d,date_min_a,date_max_a,Dvmax,ibody,accuracy)
Date0=date2mjd2000(date_min_d);%Date in mjd2000 where we start to analyze the mission
min_t_d=Date0*24*3600;%FIRST considered Departure Time in SECONDS from 1st January 2000, 12:00
max_t_d=date2mjd2000(date_max_d)*24*3600;%LAST considered Departure Time in SECONDS from 1st January 2000, 12:00
min_t_a=date2mjd2000(date_min_a)*24*3600;%FIRST considered Arrival Time in SECONDS from 1st January 2000, 12:00
max_t_a=date2mjd2000(date_max_a)*24*3600;%LAST considered Arrival Time in SECONDS from 1st January 2000, 12:00
t_d=linspace(0,max_t_d-min_t_d,accuracy);%Vector of Considered Times of DEPARTURE in SECONDS, Starting at 0 secs
t_a=linspace(min_t_a-min_t_d,max_t_a-min_t_d,accuracy);%Vector of Considered Times of ARRIVAL in SECONDS, Relative to times of departure
[kep1,mi] = uplanet(Date0, 3);
[kep2,~] = uplanet(Date0, ibody);
orbital_parameters_1.a=kep1(1);
orbital_parameters_1.ecc=kep1(2);
orbital_parameters_1.INCLI=kep1(3);
orbital_parameters_1.RAAN=kep1(4);
orbital_parameters_1.PA=kep1(5);
orbital_parameters_1.theta=kep1(6);
orbital_parameters_2.a=kep2(1);
orbital_parameters_2.ecc=kep2(2);
orbital_parameters_2.INCLI=kep2(3);
orbital_parameters_2.RAAN=kep2(4);
orbital_parameters_2.PA=kep2(5);
orbital_parameters_2.theta=kep2(6);

[Dvtot,~,TOF]=porkchop2(t_d,t_a,mi,orbital_parameters_1,orbital_parameters_2,Dvmax);
Tof=TOF/24/3600;
T_d=t_d./(24*3600);
T_a=t_a./(24*3600);
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
%and finally, we create the vector needed as the contour doc says, rounding the days to the unity.
vt=round(linspace(startTOF,endTOF,n));
figure(1);surface(T_d,T_a,Dvtot');
figure(2);contour(T_d,T_a,Dvtot',50);
hold on;contour(T_d,T_a,Tof',vt,'LineColor','k','ShowText','on');
[Dvtotmin,I] = min(Dvtot(:));
[i,j] = ind2sub(size(Dvtot),I);
T_min_d=T_d(i);
T_min_a=T_a(j);

