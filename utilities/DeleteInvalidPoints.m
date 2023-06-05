function pointArray = DeleteInvalidPoints(pointArray)
I = [];
for i = 1:size(pointArray,1)
     %% 无效点为（2, 2），其他有效状态点，因为已经归一化到[0,1]之间，保证不会被删掉
%     if (abs(pointArray(i,1)-2)<1e-5 && abs(pointArray(i,2)-2)<1e-5)  

    %% 无效点为坐标小于0
    if (pointArray(i,1)<0 || pointArray(i,2)<0) 
        I(end+1) = i;
    end
end
pointArray(I,:)=[];

end