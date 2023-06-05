function realAction = AntiNormalize(normAction, RANGE_polygon)

ring_xmin = RANGE_polygon(1);
ring_xmax = RANGE_polygon(2);
ref_x = ring_xmax - ring_xmin;

ring_ymin = RANGE_polygon(3);
ring_ymax = RANGE_polygon(4);
ref_y = ring_ymax - ring_ymin;

ref_d = max([ref_x, ref_y]);

realAction(1) = normAction(1) * ref_d + ring_xmin;
realAction(2) = normAction(2) * ref_d + ring_ymin;

end

