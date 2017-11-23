function F = hyp_system(x,v_inf_plus,v_inf_minus,delta,mi_earth)

F(1) = x(1)-1-(x(5)*(v_inf_minus)^2)/mi_earth;
F(2) = x(2)-1-(x(5)*(v_inf_plus)^2)/mi_earth;
F(3) = x(3)-2*asin(1/x(1));
F(4) = x(4)-2*asin(1/x(2));
F(5) = delta-x(3)/2-x(4)/2;

end