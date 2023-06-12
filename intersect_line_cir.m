function [flag_intersect_v] = intersect_line_cir(P1,P2,C,D_f,antigen_r)
radius = (D_f/2)+antigen_r;
radii = zeros(size(C,1),1)+radius;
no_collagen = isempty(C);
flag_intersect_v = zeros(size(C,1),1);
dis_thresh = D_f;
for i = 1:size(C,1)
    if (norm(C(i,:)-P2)<=dis_thresh)
        d = P2' - P1';
        d_uni = d/(1E+2);
        Point_array = [P1(1):d_uni(1):P2(1) ; P1(2):d_uni(2):P2(2)];
        dist_summ = ((Point_array(1,:)-C(i,1)).^2+(Point_array(2,:)-C(i,2)).^2).^0.5;

        if (sum(dist_summ<=(radii(i)+antigen_r)) ~= 0)
            flag_intersect = 1;
        else
            flag_intersect = 0;
        end

    else
        flag_intersect = 0;
    end
            flag_intersect_v(i) = flag_intersect;
end
    if no_collagen
        flag_intersect_v(i) = 0;
    end
end

