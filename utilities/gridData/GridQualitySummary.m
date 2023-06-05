function [AreaRatioQuality, ShapeQuality, AngleQuality, Skewness] = GridQualitySummary(cellNodeTopology, xCoord, yCoord, Grid_stack, faces)
if nargin == 3
    Grid_stack = ConstructGridStack(cellNodeTopology,xCoord, yCoord);
    faces = 1:size(Grid_stack,1);
end

if nargin == 4
    faces = 1:size(Grid_stack,1);
end

if isempty(cellNodeTopology)
%     cellNodeTopology = ConstructGridTopo(Grid_stack);
    cellNodeTopology = ConstructCellTopology(Grid_stack);
    PLOT(Grid_stack, xCoord, yCoord);
end

nFaces = size(faces,2);
nCells = size(cellNodeTopology,1);

if(size(cellNodeTopology,2) == 3)
    cellNodeTopology(:,4) = -1;
end
%% 计算网格角度质量
[AngleQuality, Skewness, ShapeQuality] = GridAngleQualitySummary(cellNodeTopology, xCoord, yCoord);

%% 计算Area Ratio
% AreaRatioQuality = GridAreaRatioQualitySummary(Grid_stack, faces, cellNodeTopology, xCoord, yCoord);
AreaRatioQuality = [0 0 0];
end
