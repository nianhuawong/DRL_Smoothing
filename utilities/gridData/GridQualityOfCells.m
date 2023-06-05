function [ShapeQuality, AngleQuality, Skewness] = GridQualityOfCells(cellContainer, xCoord, yCoord)
nCells = size(cellContainer,1);
minQ = 10000;
maxQ = 0;
averQ = 0;
for i = 1:nCells
    cell = cellContainer(i,:);
    if(length(cell)==3)
        quality = QualityCheckQuad(cell(1), cell(2), cell(3), -1, xCoord, yCoord);
    elseif(length(cell)==4)
        quality = QualityCheckQuad(cell(1), cell(2), cell(3), cell(4), xCoord, yCoord);
    else
        error("GridQualityOfCells，单元类型不对！");
    end
    
    minQ = min([minQ, quality]);
    maxQ = max([maxQ, quality]);
    averQ = averQ + quality;
end

averQ = averQ / nCells;

ShapeQuality = [minQ; maxQ; averQ];
[AngleQuality, Skewness] = GridAngleQualitySummary(cellContainer, xCoord, yCoord);

end