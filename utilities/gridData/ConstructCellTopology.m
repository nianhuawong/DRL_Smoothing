function cell2nodeTopo = ConstructCellTopology(Grid_stack)
nFaces = size(Grid_stack,1);
cell2node = cell(100000,1);
for i= 1: nFaces
    node1 = Grid_stack(i,1);
    node2 = Grid_stack(i,2);
    leftCell = Grid_stack(i,3);
    rightCell = Grid_stack(i,4);
    
    cell2node{leftCell} = [cell2node{leftCell}, node1, node2];
    if rightCell>0
        cell2node{rightCell} = [cell2node{rightCell}, node1, node2];
    end
end

cellCount =0;
for i=1:size(cell2node,1)
    if isempty(cell2node{i})
        continue;
    end
    
    cellCount = cellCount + 1;
    
    cell2node{i} = unique(cell2node{i}, 'stable');
    
    if(length(cell2node{i})==4)
        node1 = cell2node{i}(1);
        node2 = cell2node{i}(2);
        node3 = cell2node{i}(3);
        node4 = cell2node{i}(4);
        neighbor1 = NeighborNodes(node1, Grid_stack, -1);
        if sum(neighbor1==node3) ~= 0
            cell2node{i} = [node1, node2, node4, node3];
        elseif sum(neighbor1==node4) ~= 0
            cell2node{i} = [node1, node2, node3, node4];
        end
    end
end

cell2nodeTopo = zeros(cellCount,4);
for i=1:cellCount
    cell2nodeTopo(i,1:3) = cell2node{i}(1:3);
    cell2nodeTopo(i,4) = 0;
    if length(cell2node{i}) == 4
        cell2nodeTopo(i,4) = cell2node{i}(4);
    end
end
end