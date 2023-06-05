function InitializeRingFacesAndNeighbors(this)
%             this.ringFacesOfEachNode = zeros(this.nNodes,50);
%             this.neighborNodesOfEachNode = zeros(this.nNodes, this.max_surrounding);
this.ringFacesOfEachNode = [];
this.neighborNodesOfEachNode = [];
this.neighborCellsOfEachNode = [];
for iNode = 1:this.nNodes
    neighbors = NeighborNodes(iNode, this.Grid_stack, -1);
    [neighbors,faces] = SortStatePoints(this, neighbors);
    neighCells = ComputeNeighborCells(this, iNode);
    
    this.ringFacesOfEachNode(end+1,1:size(faces,2)) = faces;
    this.neighborNodesOfEachNode(end+1,1:size(neighbors,2)) = neighbors;
    this.neighborCellsOfEachNode(end+1,1:size(neighCells,2)) = neighCells;
end
end