function [tof,v_transf_1,v_transf_2] = pork_chopLEG2(X1,X2,tdep,tarr,misun)
l=length(tdep);
m=length(tarr);
v_transf_1 = zeros(l,m,3);
v_transf_2 = zeros(l,m,3);
tof = zeros(l,m);
for i = 1:l
    for j = 1:m
        tof(i,j) = tarr(j)-tdep(i);
        if tof(i,j) <= 0
            tof(i,j) = NaN;
            for k=1:3
                v_transf_1(i,j,k)=NaN;
                v_transf_2(i,j,k)=NaN;
            end
        else
            [~,~,~,~,v_transf_1(i,j,:),v_transf_2(i,j,:),t_p,~] = lambertMR(X1(i,1:3),X2(j,1:3),tof(i,j),misun,0,0,0,0);
            if tof(i,j) <= t_p
                tof(i,j) = NaN;
                for k=1:3
                    v_transf_1(i,j,k)=NaN;
                    v_transf_2(i,j,k)=NaN;
                end
            end
        end
    end
end
