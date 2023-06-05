classdef Segment
    properties
        SP;
        EP;
		direction;
		normal;
        len;
    end
    
    methods
         function this = Segment(P1, P2)
             this.SP = P1;
             this.EP = P2;
             
             this.direction = this.EP - this.SP;
             this.len = norm(this.direction);
             this.direction = this.direction / this.len;

             if length(P1) == 2 && length(P2) == 2
                 this.normal = [this.direction(2), -this.direction(1)];
             end
         end
         
         function [dist, tarPoint] = NearestPointOnSegment(this, point)
             AP = point - this.SP;
             projValue1 = dot(AP, this.direction);
             
             if projValue1 <= 0
                 tarPoint = this.SP;
                 dist = sqrt(dot(AP, AP));
                 return;
             end

             BP = point - this.EP;
             projValue2 = dot(BP, this.direction);
             
             if projValue2 <= 0
                 tarPoint = this.EP;
                 dist = sqrt(dot(BP, BP));
                 return;
             end             
             
             tarPoint = this.SP + projValue1 * this.direction;
             
             PQ = tarPoint - point;
             dist = sqrt(dot(PQ, PQ));
         end
         
         function [flag, point] = CrossCheckWithSegment(this, seg)
%              eps = -1e-15;
             eps = 0;
             
             a = this.SP;
             b = this.EP;
             c = seg.SP;
             d = seg.EP;
             
             % 1.计算cd的法线向量n2
             cd = d - c;
             n2 = [cd(2), -cd(1)];
             
             % 2.计算点c、点a、点b在法向量n2上的投影距离
             dist_a_n2 = dot(a, n2);
             dist_c_n2 = dot(c, n2);
             dist_b_n2 = dot(b, n2);
             
             % 3.判断是否不相交，如果a, b都在c的同一侧，则不相交
             if (dist_a_n2-dist_c_n2) * (dist_b_n2-dist_c_n2) > eps
                 flag = false;
                 point = [nan, nan];
                 return;
             end
             
             % 4.计算ab的法线向量n1
             ab = b - a;
             n1 = [ab(2), -ab(1)];
             
             % 5.计算点a、点c、点d在法向量n1上的投影距离
             dist_a_n1 = dot(a, n1);
             dist_c_n1 = dot(c, n1);
             dist_d_n1 = dot(d, n1);
             
             % 6.判断是否不相交，如果c, d都在a的同一侧，则不相交
             if (dist_c_n1-dist_a_n1) * (dist_d_n1-dist_a_n1) >eps
                 flag = false;
                 point = [nan, nan];
                 return;
             end
             
             % 7.求交点
             n1(3) = 0;
             n2(3) = 0;
             denominator = cross(n1, n2);
             denominator = denominator(3);
             fraction = (dist_c_n2 - dist_a_n2) / denominator;
             
             point = a + fraction * ab;
             flag = true;
         end
         
    end
end