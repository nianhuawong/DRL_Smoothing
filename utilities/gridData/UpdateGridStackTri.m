function Grid_stack = UpdateGridStackTri(Grid_stack, node1, node2, node3, nCells , xCoord, yCoord)

dist12 = DISTANCE(node1, node2, xCoord, yCoord);
dist13 = DISTANCE(node1, node3, xCoord, yCoord);
dist23 = DISTANCE(node2, node3, xCoord, yCoord);
%%
flag1 = IsLeftCell(node1, node2, node3, xCoord, yCoord);
[direction, row] = FrontExist(node1,node2, Grid_stack);    
if( flag1 == 1 )        %如果为左单元        
    if( row ~= -1 )  %如果已经存在
        if(direction == 1)  %如果AFT_stack_sorted中存储方向为（node1，node2）
            Grid_stack(row, 3) = nCells;%则更新到左单元
        elseif( direction == -1 )
            Grid_stack(row, 4) = nCells; %否则更新到右单元
        end                                                   
    else %如果不存在，则按反的逻辑新增阵面，即阵面方向反过来（node2, node1），nCells_AFT更新到右单元，左单元不更新
        Grid_stack(end+1,:) = [node2, node1, 0, nCells, dist12, size(Grid_stack,1)+1, 2];                     
    end
else  %如果为右单元
    if( row ~= -1 ) %如果已经存在
        if(direction == 1) %如果AFT_stack_sorted中存储方向为（node1，node2）
            Grid_stack(row, 4) = nCells; %则更新到右单元
        elseif(direction == -1)
            Grid_stack(row, 3) = nCells; %否则更新到左单元
        end                        
    else%如果不存在，则按反的逻辑新增阵面，即阵面方向反过来（node2, node1），nCells_AFT更新到左单元，右单元不更新
        Grid_stack(end+1,:) = [node2, node1, nCells, 0, dist12, size(Grid_stack,1)+1, 2];                     
    end
end

flag2 = IsLeftCell(node2, node3, node1, xCoord, yCoord);
[direction, row] = FrontExist(node2,node3, Grid_stack); 
if( flag2 == 1 )
    if( row ~= -1 )
        if(direction == 1)
            Grid_stack(row, 3) = nCells;
        elseif( direction == -1 )
            Grid_stack(row, 4) = nCells;
        end  
    else
        Grid_stack(end+1,:) = [node3, node2, 0, nCells, dist23, size(Grid_stack,1)+1, 2];                     
    end                
else
    if( row ~= -1 )
        if(direction == 1)
            Grid_stack(row, 4) = nCells; 
        elseif(direction == -1)
            Grid_stack(row, 3) = nCells;
        end 
    else
        Grid_stack(end+1,:) = [node3, node2, nCells, 0, dist23, size(Grid_stack,1)+1, 2];                     
    end
end

flag3 = IsLeftCell(node3, node1, node2, xCoord, yCoord);
[direction, row] = FrontExist(node3,node1, Grid_stack); 
if( flag3 == 1 )
    if( row ~= -1 )
        if(direction == 1)
            Grid_stack(row, 3) = nCells;
        elseif( direction == -1 )
            Grid_stack(row, 4) = nCells;
        end  
    else
        Grid_stack(end+1,:) = [node1, node3, 0, nCells, dist13, size(Grid_stack,1)+1, 2];                     
    end                
else
    if( row ~= -1 )
        if(direction == 1)
            Grid_stack(row, 4) = nCells; 
        elseif(direction == -1)
            Grid_stack(row, 3) = nCells;
        end 
    else
        Grid_stack(end+1,:) = [node1, node3, nCells, 0, dist13, size(Grid_stack,1)+1, 2];                     
    end
end    