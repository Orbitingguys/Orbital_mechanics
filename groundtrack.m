function [alpha,delta,lambda,phi,index] = groundtrack(R_geo,omega_planet,T)

N = length(R_geo(:,1));
delta = zeros(1,N);
alpha = zeros(1,N);
lambda = zeros(1,N);

for i = 1:N
    norm_r = norm(R_geo(i,:));
    delta(i) = asin(R_geo(i,3)/norm_r);
    
    if R_geo(i,2)>=0
        alpha(i) = acos(R_geo(i,1)/(norm_r*cos(delta(i))));
    else
        alpha(i) = 2*pi - acos(R_geo(i,1)/(norm_r*cos(delta(i))));
    end
    
    lambda(i) = alpha(i) - omega_planet * (T(i)-T(1));
    
    if lambda(i) < -pi
        lambda(i) = lambda(i) + 2*pi;
    elseif lambda(i) > pi
        lambda(i) = lambda(i) - 2*pi;
    end
    
           
end


delta = delta.*(180/pi);
alpha = alpha.*(180/pi);
lambda = lambda.*(180/pi);
phi = delta;

index = [];
for i = 2:N
    if (lambda(i)-lambda(i-1))< 0
        index = [index,i-1];
    end
end
% index = sort(index);

% for i = 1:length(index)
%     if index(i) == 0
%        
%        index(i) = [];
%     end
% end
    
    
    
    
    
    
    
