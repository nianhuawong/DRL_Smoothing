function normedPoint = Normalize(realPoint, RANGE_polygon)
ring_xmin = RANGE_polygon(1);
ring_xmax = RANGE_polygon(2);
ref_x = ring_xmax - ring_xmin;

ring_ymin = RANGE_polygon(3);
ring_ymax = RANGE_polygon(4);
ref_y = ring_ymax - ring_ymin;

ref_d = max([ref_x, ref_y]);

normedPoint(1) = (realPoint(1) - ring_xmin) / ref_d;
normedPoint(2) = (realPoint(2) - ring_ymin) / ref_d;
end