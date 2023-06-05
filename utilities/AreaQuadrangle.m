function area = AreaQuadrangle(node1, node2, node3, node4)
area1 = AreaTriangle(node1, node2, node3);
area2 = AreaTriangle(node1, node3, node4);
area = area1 + area2;
end