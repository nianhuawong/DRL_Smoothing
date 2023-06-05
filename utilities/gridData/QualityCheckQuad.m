function quality = QualityCheckQuad(node1, node2, node3, node4, xCoord, yCoord)
PA = [xCoord(node1), yCoord(node1)];
PB = [xCoord(node2), yCoord(node2)];
PC = [xCoord(node3), yCoord(node3)];

if node4 <= 0 
    quality = QualityCheckTri(PA, PB, PC);
    return;
end

flagConvexPoly = IsConvexPloygon(node1, node2, node3, node4, xCoord, yCoord);
if ~flagConvexPoly
    quality = 0;
    return;
end

PD = [xCoord(node4), yCoord(node4)];
%% 求对角线AC和BD的交点O
seg1 = Segment(PA, PC);
seg2 = Segment(PB, PD);

[flag, PO] = seg1.CrossCheckWithSegment(seg2);
if ~flag
    PA,PB,PC,PD
    disp("对角线AC和BD不相交");
    quality = 0;
    return;
end

%%
alpha1 = QualityCheckTri(PA, PB, PO);
alpha2 = QualityCheckTri(PB, PC, PO);
alpha3 = QualityCheckTri(PC, PD, PO);
alpha4 = QualityCheckTri(PD, PA, PO);
tmp = [alpha1,alpha2,alpha3,alpha4];
tmp = sort(tmp, 'ascend');
quality = tmp(1) * tmp(2) / tmp(3) / tmp(4);
end

function quality = QualityCheckTri(node1, node2, node3)
a = DISTANCE2(node1, node2);
b = DISTANCE2(node2, node3);
c = DISTANCE2(node3, node1);   
area = AreaTriangle(node1, node2, node3);             %三角形面积
%% 
% r = 2.0 * area / ( ( a + b + c ) );                 %内切圆半径
% R = a * b * c / 4.0 / area;                         %外接圆半径  
% quality =   3.0 * r / R;

%% 三角形网格质量的另一种求法
quality = 4.0 * sqrt(3.0) * area / ( a * a + b * b + c * c );

end

function dist = DISTANCE2(node1, node2)
vec = node1 - node2;
dist = sqrt(dot(vec, vec));
end

%% 另一种确定四边形网格质量的方法，如下
% % syms x y;
% x1 = xCoord(node1);y1 = yCoord(node1);
% x2 = xCoord(node2);y2 = yCoord(node2);
% x3 = xCoord(node3);y3 = yCoord(node3);
% x4 = xCoord(node4);y4 = yCoord(node4);
% % % % solve((y-y1)/(y3-y1)-(x-x1)/(x3-x1),(y-y2)/(y4-y2)-(x-x2)/(x4-x2), x, y);
% x = (x1*x2*y3 - x2*x3*y1 - x1*x2*y4 + x1*x4*y2 - x1*x4*y3 + x3*x4*y1 + x2*x3*y4 - x3*x4*y2)/(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3);
% y = (x1*y2*y3 - x3*y1*y2 - x2*y1*y4 + x4*y1*y2 - x1*y3*y4 + x3*y1*y4 + x2*y3*y4 - x4*y2*y3)/(x1*y2 - x2*y1 - x1*y4 + x2*y3 - x3*y2 + x4*y1 + x3*y4 - x4*y3);
% 
% node5 = length(xCoord) + 1;
% xCoord_tmp = [xCoord;x];
% yCoord_tmp = [yCoord;y];