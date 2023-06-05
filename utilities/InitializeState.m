function [initialState, RANGE] = InitializeState(ringNodesCoord)
ring_xmin = min(ringNodesCoord(:,1));
ring_xmax = max(ringNodesCoord(:,1));
ref_x = ring_xmax - ring_xmin;

ring_ymin = min(ringNodesCoord(:,2));
ring_ymax = max(ringNodesCoord(:,2));
ref_y = ring_ymax - ring_ymin;

ref_d = max([ref_x, ref_y]);

ringNodesCoord(:,1) = ringNodesCoord(:,1) - ring_xmin;
ringNodesCoord(:,2) = ringNodesCoord(:,2) - ring_ymin;

initialState = ringNodesCoord ./ ref_d;

RANGE = [ring_xmin, ring_xmax, ring_ymin, ring_ymax];
end

