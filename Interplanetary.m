function[dVtot_opt,tof_tot_opt,tau_a,tau_d]=Interplanetary(date_dep_min,date_dep_max,date_arr_min,date_arr_max,body_dep,body_fb,body_arr,I,J,K)
%% Constant definition

misun=astroConstants(4);
midep=astroConstants(body_dep+10);
miarr=astroConstants(body_arr+10);
mifb=astroConstants(body_fb+10);
RPfbmin=70000;

%% Times definition

mjd2000_dep_min_abs=date2mjd2000(date_dep_min);
mjd2000_dep_max_abs=date2mjd2000(date_dep_max);
mjd2000_arr_min_abs=date2mjd2000(date_arr_min);
mjd2000_arr_max_abs=date2mjd2000(date_arr_max);
mjd2000_dep_abs=linspace(mjd2000_dep_min_abs,mjd2000_dep_max_abs,I);
mjd2000_arr_abs=linspace(mjd2000_arr_min_abs,mjd2000_arr_max_abs,J);
mjd2000_fb_abs=linspace(mjd2000_dep_min_abs,mjd2000_arr_max_abs,K);
tdep_abs=mjd2000_dep_abs.*(24*3600);
tarr_abs=mjd2000_arr_abs.*(24*3600);
tfb_abs=mjd2000_fb_abs.*(24*3600);
t0=tdep_abs(1);
tdep=tdep_abs-t0;
tarr=tarr_abs-t0;
tfb=tfb_abs-t0;

%% Initial conditions

[orbital_parameters_dep0]=OrbitalParameters(mjd2000_dep_abs(1),body_dep);
[orbital_parameters_fb0]=OrbitalParameters(mjd2000_dep_abs(1),body_fb);
[orbital_parameters_arr0]=OrbitalParameters(mjd2000_dep_abs(1),body_arr);
[R_dep0,V_dep0] = kep2geo (orbital_parameters_dep0,misun);
[R_arr0,V_arr0] = kep2geo (orbital_parameters_arr0,misun);
[R_fb0,V_fb0] = kep2geo (orbital_parameters_fb0,misun);
Xdep0=[R_dep0;V_dep0];
Xarr0=[R_arr0;V_arr0];
Xfb0=[R_fb0;V_fb0];

%% Possition of planets

options = odeset('Reltol',1e-13,'Abstol',1e-14);

[~,Xdep] = ode113(@orbit_dynamics,tdep,Xdep0,options,misun);
Rdep=(Xdep(:,1:3))';
Vdep=(Xdep(:,4:6))';
[~,Xarr] = ode113(@orbit_dynamics,tarr,Xarr0,options,misun);
Rarr=(Xarr(:,1:3))';
Varr=(Xarr(:,4:6))';
[~,Xfb] = ode113(@orbit_dynamics,tfb,Xfb0,options,misun);
Rfb=(Xfb(:,1:3))';
Vfb=(Xfb(:,4:6))';

%% First Leg

[tof1,V_start1,V_end1] = pork_chopLEG2(Xdep,Xfb,tdep,tfb,misun);
dVdep=zeros(I,K);
for i=1:I
    for k=1:K
        dVdep(i,k)=norm(permute(V_start1(i,k,:),[3,1,2])-Vdep(:,i));
    end
end

%% Second Leg

[tof2,V_start2,V_end2] = pork_chopLEG2(Xfb,Xarr,tfb,tarr,misun);
dVarr=zeros(K,J);
for k=1:K
    for j=1:J
        dVarr(k,j)=norm(Varr(:,j)-permute(V_end2(k,j,:),[3,1,2]));
    end
end

%% Fly-by
for k=1:K
    [orbital_parameters_planet]=car2kep(Rfb(:,k),Vfb(:,k),misun);
    for i=1:I
        for j=1:J
            V_i=permute(V_end1(i,k,:),[3,1,2]);
            V_f=permute(V_start2(k,j,:),[3,1,2]);
            if any(isnan(V_i))==0 && any(isnan(V_f))==0
                [~, ~, dVfb(i,j,k)]= poweredflyby (V_i,V_f,orbital_parameters_planet,misun,mifb,RPfbmin);
                dVtot(i,j,k)=dVdep(i,k)+dVarr(k,j)+dVfb(i,j,k);
                tof_tot(i,j,k)=tof1(i,k)+tof2(k,j);
            else
                dVtot(i,j,k)=NaN;
                tof_tot(i,j,k)=tof1(i,k)+tof2(k,j);
            end
        end
    end
end
%% Optimal choice

dVtot_opt=zeros(I,J);
Kopt=zeros(I,J);
tof_tot_opt=zeros(I,J);
for i=1:I
    for j=1:J
        [dVtot_opt(i,j),Kopt(i,j)]=min(dVtot(i,j,:));
        tof_tot_opt(i,j)=tof_tot(i,j,Kopt(i,j));
    end
end
%% PlOTS

%We have to define the dates of the axes,
%We first define the tau for each time of both vectors

for i=1:length(mjd2000_dep_abs)
    tau_d(i)=mjd20002tau(mjd2000_dep_abs(i));
end
for j=1:length(mjd2000_arr_abs)
    tau_a(j)=mjd20002tau(mjd2000_arr_abs(j));
end

% 3D surface plot
figure                                      
surf(tau_d,tau_a,dVtot_opt','EdgeColor','none');
%Settings for surface 3D surface plot
title('Pork Chop Surface','FontSize',13)
xlabel('Departure time')
ylabel('Arrival time')
zlabel('Delta V optimalized')
tickDeparture=linspace(tau_d(1),tau_d(end),20);
tickArrival=linspace(tau_a(1),tau_a(end),20);
set(gca,'xtick',tickDeparture,'XTickLabelRotation',45);
set(gca,'ytick',tickArrival,'YTickLabelRotation',45);
datetick('x',1,'keepticks');
datetick('y',1,'keepticks');

% contour plot

[vt]=valuesTOF(mjd2000_dep_abs,mjd2000_arr_abs);
mjd2000of=tof_tot_opt./(3600*24);
figure                                      
contour(tau_d,tau_a,dVtot_opt',I);
hold on
contour(tau_d,tau_a,mjd2000of',vt,'LineColor','k','ShowText','on');
%Settings for contour plot
colorbar;
caxis([0,40]);
title('Pork Chop Contour','FontSize',13)
xlabel('Departure time')
ylabel('Arrival time')
set(gca,'xtick',tickDeparture,'XTickLabelRotation',45);
set(gca,'ytick',tickArrival,'YTickLabelRotation',45);
datetick('x',1,'keepticks');
datetick('y',1,'keepticks');
colorbar;
