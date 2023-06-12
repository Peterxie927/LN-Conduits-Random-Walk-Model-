function [x_temp y_temp V_temp Point] = conduit_collision(x_temp,y_temp,x_p,y_p,D,antigen_r,r,vel_j)
%CONDUIT_COLLISION Summary of this function goes here
coefficients = polyfit([x_temp, x_p], [y_temp , y_p], 1);  % line passing old point and new point which is outside circle
slope = coefficients (1);
intercpt = coefficients (2);
[xout,yout] = linecirc(slope,intercpt,0,0,(D/2)-antigen_r); % find intrsection of the line with circle
d1=sqrt((xout(1)-x_temp)^2+(yout(1)-y_temp)^2); %There are two intersections so calculate the distance from
d2=sqrt((xout(2)-x_temp)^2+(yout(2)-y_temp)^2); %new point to these intersection and choose the closer one
if d1<d2
    Point(1)=xout(1);
    Point(2)=yout(1);
else
    Point(1)=xout(2);
    Point(2)=yout(2);
end
a = Point-[x_p y_p]; %Find the vector from the intersection point to starting point
n_vec = [0,0]-Point; n_hat = n_vec./norm(n_vec);
r_vec = a-2*(dot(n_hat,a))*n_hat; r_hat = r_vec./norm(r_vec);
d_a = sqrt(a(1)^2+a(2)^2);
d_b = r-d_a;

if (d_b*r_hat(1)>=abs(r)) || (d_b*r_hat(2)>=abs(r))
    distcheck = 1;
    d_b = abs(r);
else
    distcheck = 0;
    d_b = d_b;
end
V_temp = r_hat*d_b;
x_temp = Point(1)+(d_b*r_hat(1));
y_temp = Point(2)+(d_b*r_hat(2));
end

