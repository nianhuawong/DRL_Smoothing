function Grid_stack = UpdateGridStackTri(Grid_stack, node1, node2, node3, nCells , xCoord, yCoord)

dist12 = DISTANCE(node1, node2, xCoord, yCoord);
dist13 = DISTANCE(node1, node3, xCoord, yCoord);
dist23 = DISTANCE(node2, node3, xCoord, yCoord);
%%
flag1 = IsLeftCell(node1, node2, node3, xCoord, yCoord);
[direction, row] = FrontExist(node1,node2, Grid_stack);    
if( flag1 == 1 )        %���Ϊ��Ԫ        
    if( row ~= -1 )  %����Ѿ�����
        if(direction == 1)  %���AFT_stack_sorted�д洢����Ϊ��node1��node2��
            Grid_stack(row, 3) = nCells;%����µ���Ԫ
        elseif( direction == -1 )
            Grid_stack(row, 4) = nCells; %������µ��ҵ�Ԫ
        end                                                   
    else %��������ڣ��򰴷����߼��������棬�����淽�򷴹�����node2, node1����nCells_AFT���µ��ҵ�Ԫ����Ԫ������
        Grid_stack(end+1,:) = [node2, node1, 0, nCells, dist12, size(Grid_stack,1)+1, 2];                     
    end
else  %���Ϊ�ҵ�Ԫ
    if( row ~= -1 ) %����Ѿ�����
        if(direction == 1) %���AFT_stack_sorted�д洢����Ϊ��node1��node2��
            Grid_stack(row, 4) = nCells; %����µ��ҵ�Ԫ
        elseif(direction == -1)
            Grid_stack(row, 3) = nCells; %������µ���Ԫ
        end                        
    else%��������ڣ��򰴷����߼��������棬�����淽�򷴹�����node2, node1����nCells_AFT���µ���Ԫ���ҵ�Ԫ������
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