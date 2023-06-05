function Coord_opt = SpringOptimize_delaunay(triMesh,invalidCells,wallNodes,times, xCoord_old, yCoord_old)
if nargin == 4
    xCoord_old = triMesh.Points(:,1);
    yCoord_old = triMesh.Points(:,2);
end

nNodes = length(xCoord_old);
xCoord = zeros(nNodes,1);
yCoord = zeros(nNodes,1);

ti = vertexAttachments(triMesh); %	连接到顶点的三角形或四面体
tri = triMesh.ConnectivityList;

% figure;
for NN = 1:times
    for i = 1:nNodes
        if sum(wallNodes == i) ~= 0 
            xCoord(i) = xCoord_old(i);
            yCoord(i) = yCoord_old(i);            
        else          
            attachedTri = ti(i,:);
            node_id = tri(attachedTri{:},:);
            node_id = unique(node_id(:),'stable');
            nn = length(node_id);
            xCoord_tmp = 0; yCoord_tmp = 0;
            for j =1:nn
                node = node_id(j);
                xCoord_tmp = xCoord_old(node) + xCoord_tmp;
                yCoord_tmp = yCoord_old(node) + yCoord_tmp;
            end
            xCoord(i) = xCoord_tmp / nn;
            yCoord(i) = yCoord_tmp / nn;
        end
        xCoord_old(i) = xCoord(i);%更新点i坐标后，立马采用新坐标进行后续更新
        yCoord_old(i) = yCoord(i);
    end

    tri_plot = tri;
    tri_plot(invalidCells,:) = [];
    clf
    triplot(tri_plot,xCoord,yCoord);
    axis equal;
    axis off
    hold on;
    pause(0.5)
end
Coord_opt = [xCoord, yCoord];
end