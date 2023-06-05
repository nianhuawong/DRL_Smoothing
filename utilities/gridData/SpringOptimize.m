function Grid_coord = SpringOptimize(Grid_stack, BC_nodes, Grid_coord, times)

nNodes = size(Grid_coord,1);

figure;
for k=1:times
    disp(['Laplacian smooth iterations = ', num2str(k), ' / ', num2str(times)]);
    for i=1:nNodes
        
        if(sum(BC_nodes==i)>0) 
            continue;
        end
        
        neigbhors = NeighborNodes(i, Grid_stack, -1);
        
        neibghborCoords = Grid_coord(neigbhors,:);
        
        Grid_coord(i,:) = Centroid(neibghborCoords);
    end
    
    PLOT(Grid_stack, Grid_coord(:,1), Grid_coord(:,2), 0);
    pause(0.001);
end
end
