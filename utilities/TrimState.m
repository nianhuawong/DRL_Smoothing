function initialState = TrimState(initialState, max_sur)
n_sur = size(initialState,1);
if  n_sur < max_sur
    %% 无效点为坐标小于0
    initialState(n_sur+1:max_sur,:)=-1e-15*ones(max_sur-n_sur,2); 
elseif n_sur > max_sur
    initialState(max_sur+1:end, :)=[];
%     warning(['nsurroundings > ',num2str(max_sur)]);
end
end