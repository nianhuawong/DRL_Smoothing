function PLOT(AFT_stack, xCoord, yCoord, num_label)
if nargin==3
    num_label = 0;
    flag_label=[];
end
nNodes = size(xCoord,1);
flag_label = zeros(nNodes,1);

fig = clf;
% fig = figure;
fig.Color = 'white'; hold on;
len = size(AFT_stack,1);
for i = 1:len
    node1 = AFT_stack(i,1);
    node2 = AFT_stack(i,2);
    
    dist = DISTANCE(node1,node2,xCoord,yCoord);
    
    xx = [xCoord(node1),xCoord(node2)];
    yy = [yCoord(node1),yCoord(node2)];

    plot( xx, yy, '-b');  %'MarkerSize',14,'LineWidth',1
    hold on;
end

for i = 1 : nNodes
    str = num2str(i);
    if  num_label == 1 && flag_label(i) == 0
        text(xCoord(i)+0.00005*dist,yCoord(i)+0.00005*dist,str, 'Color', 'red', 'FontSize', 9)
        flag_label(i) = 1;
    end
end
axis equal;
axis off;
    
end