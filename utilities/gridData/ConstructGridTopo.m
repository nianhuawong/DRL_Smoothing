function cellNodes = ConstructGridTopo(Grid_stack)
cellNodes = [];
nCells = 100000;
cellNodes_Tmp = zeros(nCells,1);
nFaces = size(Grid_stack,1);
for i= 1: nFaces
    node1 = Grid_stack(i,1);
    node2 = Grid_stack(i,2);
    leftCell = Grid_stack(i,3);
    rightCell = Grid_stack(i,4);

    END = size(cellNodes_Tmp(leftCell,:),2) - sum(cellNodes_Tmp(leftCell,:)==0) + 1;
    cellNodes_Tmp(leftCell,END:END+1)  = [node1, node2];
    
    if rightCell ~= 0 
        END = size(cellNodes_Tmp(rightCell,:),2) - sum(cellNodes_Tmp(rightCell,:)==0) + 1;
        cellNodes_Tmp(rightCell,END:END+1) = [node1, node2];
    end
end

for j = 1:size(cellNodes_Tmp,1)
    tmp = unique(cellNodes_Tmp(j,:), 'stable');
    tmp(tmp==0) = [];
    if length(tmp) == 3
        cellNodes(j,:) = [tmp,-1];        
    elseif length(tmp) == 4
%         cellNodes(j,:) = tmp;
        node1 = tmp(1);
        node2 = tmp(2);
        node3 = tmp(3);
        node4 = tmp(4);
        neighbor1 = NeighborNodes(node1, Grid_stack, -1);
        if sum(neighbor1==node3) ~= 0 
            cellNodes(j,1) = node1;
            cellNodes(j,2) = node2;
            cellNodes(j,3) = node4;
            cellNodes(j,4) = node3;
        elseif sum(neighbor1==node4) ~= 0
            cellNodes(j,1) = node1;
            cellNodes(j,2) = node2;
            cellNodes(j,3) = node3;
            cellNodes(j,4) = node4;
        end
    else
        
    end
end
    