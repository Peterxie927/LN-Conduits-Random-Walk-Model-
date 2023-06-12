function [x_temp2 y_temp2 flag r_copy2 Point] = collagen_collision(x_temp,y_temp,antigen_prev,collagen,D_f,antigen_r,r_copy)
x_p = antigen_prev(1); y_p = antigen_prev(2);
coefficients = polyfit([x_p,x_temp], [y_p,y_temp], 1);  % line passing old point and new point which is outside circle
slope = coefficients(1);
intercpt = coefficients(2);
[xout,yout] = linecirc(slope,intercpt,collagen(1),collagen(2),(D_f/2)+antigen_r); % find intrsection of the line with circle
d1=sqrt((xout(1)-x_p)^2+(yout(1)-y_p)^2); %There are two intersections so calculate the distance from
d2=sqrt((xout(2)-x_p)^2+(yout(2)-y_p)^2); %new point to these intersection and choose the closer one
if d1<d2
    Point(1)=xout(1);
    Point(2)=yout(1);
else
    Point(1)=xout(2);
    Point(2)=yout(2);
end
a = [x_p y_p]-Point; % Find the vector from the intersection point to original vector from start of random walk
a_hat = a./norm(a);
n_vec = Point-[collagen(1),collagen(2)]; n_hat = n_vec./norm(n_vec);
r_vec = -(a_hat-2*(dot(n_hat,a_hat))*n_hat); r_hat = r_vec./norm(r_vec);
d_a = sqrt(a(1)^2+a(2)^2);
d_b = r_copy-d_a;
x_temp2 = Point(1)+(d_b*r_hat(1));
y_temp2 = Point(2)+(d_b*r_hat(2));
if (d_b*r_hat(1)>=abs(r_copy)) || (d_b*r_hat(2)>=abs(r_copy))
    distcheck = 1;
    d_b = abs(r_copy);
else
    distcheck = 0;
    d_b = d_b;
end

 x_temp2 = Point(1)+(d_b*r_hat(1));
 y_temp2 = Point(2)+(d_b*r_hat(2));
 r_copy2 = r_copy - d_a;

% V_temp = r_hat*norm(V(j,1:2));
if isnan(x_temp2) || isnan(y_temp2)
    flag = 1;
    x_temp2 = x_p;
    y_temp2 = y_p;
    r_copy2 = r_copy;
%     V = (r_copy/sqrt(3)).*randn([1,3]);
%     x_temp2 = antigen_prev(1)+V(1,1); y_temp2 = antigen_prev(2)+V(1,2); antigen(j,1:2) = [x_temp2 y_temp2];
else
    flag = 0;
end

end

