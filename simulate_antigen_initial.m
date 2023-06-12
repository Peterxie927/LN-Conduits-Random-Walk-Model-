function [antigen,c] = simulate_antigen_initial(antigen_r,N_a,D,x0,y0,collagen,N_f,D_f,L_conduit,plot_yes)
prev_x2 = zeros(N_a,1);
prev_y2 = zeros(N_a,1);
antigen = zeros(N_a,3);
c = 0; z0 = 0;
fails = 0;
prev_x = collagen(:,1);
prev_y = collagen(:,2);
for i = 1:N_a
    newAntigen = 0;
    fails = 0;
    while ~newAntigen & c<=N_a
        t = 2*pi*rand(1,1);
        r = (((D))/2)*(rand(1,1));
        x = x0 + r.*cos(t);
        y = y0 + r.*sin(t);
        z = z0 + rand(1,1)*(L_conduit-2*antigen_r)+antigen_r;
        prev_x2 = antigen(1:i-1,1);
        prev_y2 = antigen(1:i-1,2);
        prev_z2 = antigen(1:i-1,3);
        dist_summ = ((prev_x-x).^2+(prev_y-y).^2).^0.5;
        dist_summ2 = ((prev_x2-x).^2+(prev_y2-y).^2+(prev_z2-z).^2).^0.5;

        if (((x^2+y^2)^0.5)<=((D-2*antigen_r)/2))&&(((x^2+y^2)^0.5)>=((D-6*antigen_r)/2)) &&(sum(dist_summ<=((D_f/2)+(antigen_r))) == 0)
            newAntigen = 1;
            antigen(i,:) = [x y z];
            c = c+1;
            %             circle3(x,y,antigen_r,[0.6350 0.0780 0.1840]);
            fail = 0;
        else
            fails = fails + 1;
            if fails > 10000000
                disp('Cannot simulate any more antigens, Simulated antigens:');
                disp(i)
                break
            end
        end
    end
end
if plot_yes == 1
    figure (1)
    for p = 1:size(antigen,1)
        [x,y,z] = sphere;
        x1 = x*antigen_r;
        y1 = y*antigen_r;
        z1 = z*antigen_r;
        circle3(antigen(p,1),antigen(p,2),antigen_r,'b',1);
        %     plotfun = @(c,antigen_r) surf(X*antigen_r + c(1),Y*antigen_r + c(2),Z*antigen_r + c(3),...
        %     'FaceColor',.7*[1 1 1],'EdgeColor','none',...
        %     'FaceLighting','gouraud','AmbientStrength',0.5);
        %     light('Position',[-1 0 0]);
        hold on;
        grid on;
        box on
        set(gcf,'color','w')
        xlabel('X', 'FontSize', 20);
        ylabel('Y', 'FontSize', 20);
        zlabel('Z', 'FontSize', 20);
        axis equal
        hold on
    end
    xlim([-5E-07 5E-07]);
    ylim([-5E-07 5E-07]);
    zlim([0E-07 10E-07]);
end
end

