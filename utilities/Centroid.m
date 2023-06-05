function centroid = Centroid(ringNodesCoord)

x_mid = sum(ringNodesCoord(:,1))/size(ringNodesCoord,1);
y_mid = sum(ringNodesCoord(:,2))/size(ringNodesCoord,1);

centroid = [x_mid, y_mid];
end