function [AngleQuality, Skewness, ShapeQuality] = GridAngleQualitySummary(cellContainer, xCoord, yCoord)

nCells = size(cellContainer,1);
%%
minAngle = 1000;
averMinAngle = 0;
maxAngle = 0;
averMaxAngle = 0;
%%
minQuality = 1000;
maxQuality = 0;
averQuality = 0;
pos_minQ = cell(1,1);
%%
CellSkewness = zeros(nCells,1);
CellMaxAngle = zeros(nCells,1);
CellMinAngle = ones(nCells,1) * 1000;
for i=1:nCells
    CELL = cellContainer(i,:);
    CELL(CELL<=0)=[];
    
    p1 = [xCoord(CELL(1)), yCoord(CELL(1))];
    p2 = [xCoord(CELL(2)), yCoord(CELL(2))];
    p3 = [xCoord(CELL(3)), yCoord(CELL(3))];
    
    if length(CELL) == 3
        [minVal, maxVal] = SingleTriCellIncludeAngle(CELL, p1, p2, p3);
        theta_e = 60;        
        quality = QualityCheckQuad(CELL(1), CELL(2), CELL(3), -1, xCoord, yCoord);
    elseif length(CELL) == 4
        p4 = [xCoord(CELL(4)), yCoord(CELL(4))];
        [minVal, maxVal] = SingleQuadCellIncludeAngle(CELL, p1, p2, p3, p4);
        theta_e = 90;
        quality = QualityCheckQuad(CELL(1), CELL(2), CELL(3), CELL(4), xCoord, yCoord);
    end
        
    maxQuality = max([quality, maxQuality]);
    if quality < minQuality
        minQuality = quality;
        pos_minQ = i;
    end
    %minQuality = min([quality, minQuality]);
    averQuality = averQuality + quality;
    %%
    averMinAngle = averMinAngle + minVal;
    averMaxAngle = averMaxAngle + maxVal;
    
    skew1 = (theta_e - minVal)/theta_e;
    skew2 = (maxVal - theta_e)/(180 - theta_e);
    
    CellMinAngle(i) = minVal;
    CellMaxAngle(i) = maxVal;
    
    CellSkewness(i) = max([skew1, skew2]);
    
    minAngle = min([minAngle,minVal]);
    maxAngle = max([maxAngle,maxVal]);
    
    if(CellSkewness(i)>0.99)
        kkk = 1;
    end
end
averMinAngle = averMinAngle / nCells;
averMaxAngle = averMaxAngle / nCells;
averQuality = averQuality / nCells;

averMinAngle0 = mean(CellMinAngle);
averMaxAngle0 = mean(CellMaxAngle);
minAngle0 = min(CellMinAngle);
maxAngle0 = max(CellMaxAngle);

maxEquiAngleSkewness  = max(CellSkewness);
averEquiAngleSkewness = mean(CellSkewness);

AngleQuality = [minAngle;maxAngle;averMinAngle;averMaxAngle];
Skewness = [maxEquiAngleSkewness; averEquiAngleSkewness];

ShapeQuality = [minQuality; maxQuality; averQuality; pos_minQ];
end

function [minVal, maxVal] = SingleTriCellIncludeAngle(CELL, p1, p2, p3)
v12 = p2 - p1;
v13 = p3 - p1;
v12 = v12 / norm(v12);
v13 = v13 / norm(v13);
angle1 = acos(v12 * v13') * 180 / pi;

v23 = p3 - p2;
v21 = p1 - p2;
v21 = v21 / norm(v21);
v23 = v23 / norm(v23);
angle2 = acos(v23 * v21') * 180 / pi;

v31 = -v13;
v32 = -v23;
angle3 = acos(v31 * v32') * 180 / pi;

if abs((angle1+angle2+angle3)-180)>1e-4
    angle1, angle2, angle3, CELL
    error("三角形内角和不等于180°，请检查！")
end

minVal = min([angle1, angle2, angle3]);
maxVal = max([angle1, angle2, angle3]);
end

function [minVal, maxVal] = SingleQuadCellIncludeAngle(CELL, p1, p2, p3, p4)
v12 = p2 - p1;
v14 = p4 - p1;
v12 = v12 / norm(v12);
v14 = v14 / norm(v14);
angle1 = acos(v12 * v14') * 180 / pi;

v23 = p3 - p2;
v21 = p1 - p2;
v23 = v23 / norm(v23);
v21 = v21 / norm(v21);
angle2 = acos(v23 * v21') * 180 / pi;

v32 = p2 - p3;
v34 = p4 - p3;
v32 = v32 / norm(v32);
v34 = v34 / norm(v34);
angle3 = acos(v32 * v34') * 180 / pi;

v43 = p3 - p4;
v41 = p1 - p4;
v43 = v43 / norm(v43);
v41 = v41 / norm(v41);
angle4 = acos(v43 * v41') * 180 / pi;

if abs((angle1+angle2+angle3+angle4)-360)>1e-4
%     angle1, angle2, angle3, angle4, CELL
    disp("四边形内角和不等于360°，非凸四边形，请检查！")
    minVal = 0;
    maxVal = 180;    
else
    
    minVal = min([angle1, angle2, angle3, angle4]);
    maxVal = max([angle1, angle2, angle3, angle4]);
end
end