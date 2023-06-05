function [neighborCellsIndex, neighborCells] = ComputeNodeNeighborCells(cellNodeTopology, nodeIndex)
neighborCells = [];
neighborCellsIndex=[];
nCells = size(cellNodeTopology,1);
for i=1:nCells
	cell = cellNodeTopology(i,:);
    for j = 1:length(cell)
        node = cell(j);
        if(node == nodeIndex)
			neighborCells(end+1,:) = cell;
            neighborCellsIndex(end+1) = i;
            break;
        end
    end
end
end