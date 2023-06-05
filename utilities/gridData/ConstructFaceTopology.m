function Grid_stack = ConstructFaceTopology(cell2node, xCoord, yCoord)

nCells = size(cell2node,1);
Grid_stack = [];

for i= 1:nCells
    node1 = cell2node(i,1);
    node2 = cell2node(i,2);
    node3 = cell2node(i,3);
    
    Grid_stack = UpdateGridStackTri(Grid_stack, node1, node2, node3, i, xCoord, yCoord);
end
end