function flag = SinglePointOutOfDomain(newpoint, RANGE, pointArray)
if nargin <= 2
    pointArray = [];
end

% figure
% for i = 1:size(pointArray,1)
%     if i == size(pointArray,1)
%         next = 1;
%     else
%         next = i+1;
%     end
%     plot([pointArray(i,1),pointArray(next,1)],[pointArray(i,2),pointArray(next,2)],'r-o');
%     hold on;
% end
% plot(newpoint(1),newpoint(2),'g*')
% hold off;

testx = newpoint(1);
testy = newpoint(2);

Xmin = RANGE(1);
Xmax = RANGE(2);
Ymin = RANGE(3);
Ymax = RANGE(4);

flag = false;
if testx < Xmin || testx > Xmax || testy < Ymin || testy > Ymax
    flag = true;
    return;
end

vertx = pointArray(:,1);
verty = pointArray(:,2);
nCorners = size(pointArray,1);

%%
% PNPOLY - Point Inclusion in Polygon Test
% W. Randolph Franklin (WRF) https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
%% 遍历每一条边,判断是否在多边形内部 
j = nCorners;
for i=1:nCorners
    if((verty(i)>testy) ~= (verty(j)>testy)) ...
            && (testx < (vertx(j)-vertx(i)) * (testy-verty(i)) / (verty(j)-verty(i)) + vertx(i))
        
        flag = not(flag);
    end
    j = i;
end

%% 是否在外部
flag = not(flag);
end