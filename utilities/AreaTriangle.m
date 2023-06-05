function area = AreaTriangle(node1, node2, node3)
v12 = node2-node1;
v13 = node3-node1;

if length(v12) == 2 ||  length(v13) == 2
    v12(3) = 0.0;
    v13(3) = 0.0;
end

tmp = cross(v12, v13);
area = 0.5 * norm(tmp);
end