function PlotCenter(State, RANGE_polygon)
hold on;

% State = DeleteInvalidPoints(State);

nRingNodes = size(State,1);
surrounding_points = zeros(nRingNodes,2);

for i =1:nRingNodes
    surrounding_points(i,:) = AntiNormalize(State(i,:), RANGE_polygon);
end

centroid = Centroid(surrounding_points);

plot(centroid(1), centroid(2), 'ro','MarkerSize',8)
end