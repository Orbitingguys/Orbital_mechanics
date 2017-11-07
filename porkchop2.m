function [Dvtot,TPAR,TOF]=porkchop2(t_d,t_a,mi,orbital_parameters_1,orbital_parameters_2,Dvmax)
a_1 = orbital_parameters_1.a;
ecc_1 = orbital_parameters_1.ecc;
RAAN_1 = orbital_parameters_1.RAAN;
PA_1 = orbital_parameters_1.PA;
INCLI_1 = orbital_parameters_1.INCLI;
theta0_1= orbital_parameters_1.theta;
a_2 = orbital_parameters_2.a;
ecc_2 = orbital_parameters_2.ecc;
RAAN_2 = orbital_parameters_2.RAAN;
PA_2 = orbital_parameters_2.PA;
INCLI_2 = orbital_parameters_2.INCLI;
theta0_2= orbital_parameters_2.theta;
h_1=sqrt(a_1*mi*(1-ecc_1^2));
h_2=sqrt(a_2*mi*(1-ecc_2^2));
tp1=-timesinceper(theta0_1,h_1,ecc_1,mi);
tp2=-timesinceper(theta0_2,h_2,ecc_2,mi);

m=length(t_d);
n=length(t_a);
TOF=zeros(m,n);
TPAR=zeros(m,n);
Rd=zeros(m,3);
Ra=zeros(n,3);
Dvtot=zeros(m,n);
for i=1:m
    for j=1:n
        [theta_d(i)]=theta_t(t_d(i)-tp1,mi,h_1,ecc_1);
        [Rd(i,:),V1(i,:)] = kep2geo (orbital_parameters_1,mi,theta_d(i));
        [theta_a(j)]=theta_t(t_a(j)-tp2,mi,h_2,ecc_2);
        [Ra(j,:),V2(j,:)] = kep2geo (orbital_parameters_2,mi,theta_a(j));
        TOF(i,j)=t_a(j)-t_d(i);
        Dvtot(i,j)=NaN;
        if TOF(i,j)>0
            [~,~,~,ERROR(i,j),Vd,Va,TPAR(i,j),~] = lambertMR(Rd(i,:),Ra(j,:),TOF(i,j),mi,0,0,0,0);
            if TOF(i,j)>TPAR(i,j)
                Dvd=norm(V2(j,:)-Va);
                Dva=norm(Vd-V1(i,:));
                if Dvd<=Dvmax && Dva<=Dvmax
                    Dvtot(i,j)=Dvd+Dva;
                end
            end
        end 
    end
end
