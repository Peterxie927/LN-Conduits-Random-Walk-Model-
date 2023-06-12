function h = circle3(x,y,r,color,fill)
d = r*2;
px = x-r;
py = y-r;
if fill == 1
    facecolor = [0 0 1];
    h = rectangle('Position',[px py d d],'Curvature',[1,1],'EdgeColor',color,'FaceColor',facecolor);
else
    h = rectangle('Position',[px py d d],'Curvature',[1,1],'EdgeColor',color);
end
daspect([1,1,1])
end
