function [collagen,c] = simulate_collagen(D_f,N_f,D,x0,y0,L_conduit,plot_yes)
collagen = zeros(N_f ,2);
prev_x = zeros(N_f ,1);
prev_y = zeros(N_f ,1);
c = 0;
fails = 0;
for i = 1:N_f
    newCircle = 0;
    while ~newCircle 
        t = 2*pi*rand(1,1);
        r = ((D-D_f)/2)*(rand(1,1));
        x = x0 + r.*cos(t);
        y = y0 + r.*sin(t);
        prev_x = collagen(1:i-1,1);
        prev_y = collagen(1:i-1,2);
        dist_summ = ((prev_x-x).^2+(prev_y-y).^2).^0.5;

        if ((sum(dist_summ<=1.1*D_f) == 0)&&((x^2+y^2)^0.5<=((D-1.1*D_f)/2)))&&((x^2+y^2)^0.5>=(D_f))
            newCircle = 1;
            collagen(i,:) = [x y];
            c = c+1;
%             circle3(x,y,D_f/2,[0 0.4470 0.7410]);
            fail = 0;
        else
        fails = fails + 1;
            if fails > 1000000
                disp('Cannot simulate any more collagen fibres, Simulated fibres:');
                disp(i-1)
                break
            end
        end
%     hold on
    end
end
if plot_yes == 1
    figure (1)
    for i = 1:size(collagen,1)
        [Xc,Yc,Zc] = cylinder(D_f/2,100);
        h = L_conduit;
        Zc = (Zc*h);
        circle3(collagen(i,1),collagen(i,2),D_f/2,'k',0);
        hold on
    end
end

end
