function [Grid_stack, Grid_coord, cell2node, TR]  = read_STL(boundaryFile)
BC_WALL = 3;        %BC=3为物面
BC_FAR  = 9;        %BC=9为外场边界
% BC_INTERIOR = 2;    %BC=2为内部面

TR = stlread(boundaryFile);
% triplot(TR);
% axis equal;
% axis off;
% F = freeBoundary(TR);

cell2node = TR.ConnectivityList;
Grid_coord = TR.Points(:,1:2);
Grid_stack = ConstructFaceTopology(TR.ConnectivityList, Grid_coord(:,1), Grid_coord(:,2));

for i = 1:size(Grid_stack,1)
    leftCell = Grid_stack(i,3);
    rightCell = Grid_stack(i,4);
    if(leftCell>0 && rightCell==0)
        Grid_stack(i,7) = BC_WALL;
    end
    
    if(rightCell>0 && leftCell==0)
        Grid_stack(i,7) = BC_WALL;
        
        tmp = Grid_stack(i,3);
        Grid_stack(i,3) = Grid_stack(i,4);
        Grid_stack(i,4) = tmp;
    end
end

% PLOT(Grid_stack, Grid_coord(:,1), Grid_coord(:,2));
end