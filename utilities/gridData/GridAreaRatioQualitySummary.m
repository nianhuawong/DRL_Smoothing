function AreaRatioQuality = GridAreaRatioQualitySummary(Grid_stack, faces, cellNodeTopology, xCoord, yCoord)
nFaces = size(faces,2);
areaRatio = zeros(1,nFaces);
for i = 1:nFaces
    faceIndex = faces(i);
    
    leftCell = Grid_stack(faceIndex,3);
    rightCell = Grid_stack(faceIndex,4);
    if(leftCell<0 && rightCell>0)
        leftCell = rightCell;
        rightCell = -1;
    end
    node5 = cellNodeTopology{leftCell}(1);
    node6 = cellNodeTopology{leftCell}(2);
    node7 = cellNodeTopology{leftCell}(3);
    node8 = -1;
    if length(cellNodeTopology{leftCell}) >3
        node8 = cellNodeTopology{leftCell}(4);
    end
    
    if rightCell <= 0
        areaRatio(i)      = -1;
    else
        node9  = cellNodeTopology{rightCell}(1);
        node10 = cellNodeTopology{rightCell}(2);
        node11 = cellNodeTopology{rightCell}(3);
        node12 = -1;
        if length(cellNodeTopology{rightCell}) >3
            node12 = cellNodeTopology{rightCell}(4);
        end
        
        areaL = AreaQuadangleIndex(node5, node6, node7, node8, xCoord, yCoord);
        areaR = AreaQuadangleIndex(node9, node10, node11, node12, xCoord, yCoord);
        areaRatio(i) = min([areaL,areaR]) / max([areaL,areaR]);
    end
end
%% areaRatio
II = areaRatio==-1;
areaRatio( II ) = 1000;
[minAreaRatio, pos1] = min(areaRatio);

% Edge = Grid_stack(pos1,1:2);
% plot(xCoord(Edge),yCoord(Edge),'r-*');

areaRatio( II ) = [];
maxAreaRatio = max(areaRatio);
averAreaRatio = mean(areaRatio);

%% 返回最小面积比的face编号和minQuality的cell编号
pos_minAR = faces(pos1);

AreaRatioQuality = [minAreaRatio; maxAreaRatio; averAreaRatio; pos_minAR];
end

function area = AreaQuadangleIndex(nodeIn1, nodeIn2, nodeIn3, nodeIn4, xCoord, yCoord)
node1 = [xCoord(nodeIn1), yCoord(nodeIn1)];
node2 = [xCoord(nodeIn2), yCoord(nodeIn2)];
node3 = [xCoord(nodeIn3), yCoord(nodeIn3)];
if nodeIn4 > 0
    node4 = [xCoord(nodeIn4), yCoord(nodeIn4)];
    area = AreaQuadrangle(node1, node2, node3, node4);
else
    area = AreaTriangle(node1, node2, node3);
end
end

% function area = AreaTriangleIndex(nodeIn1, nodeIn2, nodeIn3, xCoord, yCoord)
% node1 = [xCoord(nodeIn1), yCoord(nodeIn1)];
% node2 = [xCoord(nodeIn2), yCoord(nodeIn2)];
% node3 = [xCoord(nodeIn3), yCoord(nodeIn3)];
%
% area = AreaTriangle(node1, node2, node3);
% end