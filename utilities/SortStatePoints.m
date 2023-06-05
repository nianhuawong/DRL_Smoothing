function  neighbors = SortStatePoints(Grid_stack, neighbors, bcFlag)
index = [];
for i = 1:size(neighbors,2)
    for j = i+1:size(neighbors,2)
        node1 = neighbors(i);
        node2 = neighbors(j);
        [~, row] = FrontExist(node1, node2, Grid_stack);
        if row > 0
            index(end+1) = row;
        end
    end
end

if(bcFlag)
    if(size(index,2)~=(size(neighbors,2)-1))
        error("边界node的ring node排序出问题！");
    end
    
    count = 1;  %对于边界node，因为没有形成闭环，所以需要从1开始
else
    if(size(index,2)~=(size(neighbors,2)))
        error("内部node的ring node排序出问题！");
    end
    
    count = 2;
end

neighbors = Grid_stack(index(1),1:2);

while(count<size(index,2))
    i = randi([2,size(index,2)]);
    node1 = Grid_stack(index(i),1);
    node2 = Grid_stack(index(i),2);
    
    if(node1 == neighbors(end) && node2~=neighbors(end-1))
        neighbors(end+1) = node2;
        count = count + 1;
    elseif(node2 == neighbors(end)&& node1~=neighbors(end-1))
        neighbors(end+1) = node1;
        count = count + 1;
    elseif(node1==neighbors(1)&& node2~=neighbors(2))
        neighbors = [node2, neighbors];
        count = count + 1;
    elseif(node2==neighbors(1)&& node1~=neighbors(2))
        neighbors = [node1, neighbors];
        count = count + 1;
    end
end
kkk = 1;
end