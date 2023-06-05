function Coord = InitNodePositions(Grid_stack, Coord, wallNodes, coeff)
rng(now); 
% rng(1);
a = -1; b = 1;
for i = 1:size(Grid_stack,1)
    node1 = Grid_stack(i,1);
    node2 = Grid_stack(i,2);
    len = Grid_stack(i,5);
    
    if sum(node1==wallNodes) == 0
        r = a + (b-a).*rand(1,1);
        Coord(node1,1) = Coord(node1,1) + coeff * r * len;
        r = a + (b-a).*rand(1,1);
        Coord(node1,2) = Coord(node1,2) + coeff * r * len;
    end
    
    if sum(node2==wallNodes) == 0
         r = a + (b-a).*rand(1,1);
        Coord(node2,1) = Coord(node2,1) + coeff * r * len;
         r = a + (b-a).*rand(1,1);
        Coord(node2,2) = Coord(node2,2) + coeff * r * len;
    end
end
end