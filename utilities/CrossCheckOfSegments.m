function [flag, point] = CrossCheckOfSegments(seg1, seg2)
[flag, point] = seg1.CrossCheckWithSegment(seg2);
end