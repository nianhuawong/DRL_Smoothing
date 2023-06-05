classdef (Abstract) DRL_Opt_Abstract < rl.env.MATLABEnvironment
    properties
        perturb_coeff   = 0.0;     % 初始网格扰动比例
        reward_coeff    = 1.0;      % minQuality在reward中的权重
        lamda           = 0.0;      % shape quality所占权重，skewness权重为1-lamda
        
        boundaryFile = './grid/training_mesh.stl';
        
        gridType = 0;       % 0-纯三角形网格，1-三角形/四边形混合网格
        
        flag_second_neighbor = false;     % 是否取邻点的邻点作为输入 
        plotOptProcess       = false;
        validation_flag      = false;    % 是否是测试验证过程
    end    
    properties
        max_ringNodes;                  % 最大临近点数量      
        Penalty         = -1;           % 统一的Penalty
        nNodes          = 0;            % 当前网格的节点数            
        nodeIndex       = 0;            % 当前遍历的nodeIndex
        
        Coord0  = [];                   % 初始网格点坐标（加了扰动之后的）     
        Coord   = [];                   % 动作完成后的坐标，实时更新
        
        Grid_stack          = [];       % 网格面拓扑    
        BC_nodes            = [];       % 边界点编号
        cellNodeTopology    = [];       % 单元拓扑，每个单元的点的编号
        Wall_nodes          = [];

        RANGE_domain    = [];       % 整个计算域的RANGE        
        RANGE_polygon   = [];       % 单个多边形的RANGE    
    end   
    properties
        State = [];
        ring_nodes4_each_node       = [];
        second_neighbor4_each_nodes = [];
        surround_nodes4_each_node   = [];
        
        ring_nodes_index = [];
        ring_nodes_coord = [];

        surround_nodes_index = [];
        surround_nodes_coord = [];
        
        neighborCells  = [];
        gridDataObj;
        actionStorage  = [];
        lapCount = 0;
    end
    properties(Access = protected)
        IsDone = false
    end
    
    methods
        function this = DRL_Opt_Abstract(ActionInfo, max_ringNodes, boundaryFile, perturb_coeff, gridType)
            ObservationInfo = rlNumericSpec([max_ringNodes 2]);
            ObservationInfo.Name = 'mesh DRL States';
            ObservationInfo.Description = 'ring nodes coords';
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            
            this.max_ringNodes = max_ringNodes;
            
            if ~isempty(boundaryFile)
                this.boundaryFile = boundaryFile;
            end
            
            if perturb_coeff >=0
                this.perturb_coeff = perturb_coeff;
            end
            
            if gridType >= 0
                this.gridType = gridType;
            end
            
            this.gridDataObj = GridDataClass(this.boundaryFile, this.gridType);
            
%             PLOT(this.gridDataObj.Grid_stack, this.gridDataObj.Grid_coord(:,1), this.gridDataObj.Grid_coord(:,2), 1);
            
            Initialize(this);
        end
        
        function Initialize(this)
            this.gridDataObj.initNodePositions(this.perturb_coeff);           
            this.gridDataObj.ConstructRingNodes();
            this.gridDataObj.ConstructNeighborNodes();
            this.gridDataObj.ConstructCell2NodeTopology();
            
            this.Coord0     = this.gridDataObj.Grid_coord;
            this.Grid_stack = this.gridDataObj.Grid_stack;
            this.BC_nodes   = this.gridDataObj.BC_nodes;
            this.Wall_nodes = this.gridDataObj.Wall_nodes;
            
            this.ring_nodes4_each_node = this.gridDataObj.Ring_nodes;
            
            if this.flag_second_neighbor
                this.surround_nodes4_each_node = this.gridDataObj.secondNeighbor;
            else
                this.surround_nodes4_each_node = this.gridDataObj.Ring_nodes;
            end
                       
            this.cellNodeTopology = this.gridDataObj.cellNodeTopology;
            
            this.nNodes  = size(this.Coord0,1);
            Xmax = max(this.Coord0(:,1));
            Ymax = max(this.Coord0(:,2));
            Xmin = min(this.Coord0(:,1));
            Ymin = min(this.Coord0(:,2));
            this.RANGE_domain = [Xmin, Xmax, Ymin, Ymax];
        end
        
        function observation = reset(this)
            this.Coord = this.Coord0;
            
            this.nodeIndex = 0;
            NodeIndexIncrement(this);
%             PLOT(this.Grid_stack, this.Coord0(:,1), this.Coord0(:,2));
            if(~this.validation_flag && this.plotOptProcess)
                PLOT(this.Grid_stack, this.Coord0(:,1), this.Coord0(:,2));
            end       
            
            ComputeStateForCurrentNode(this, this.nodeIndex);
            
            observation = this.State;
            
            this.actionStorage = zeros(this.nNodes,2);
        end
                
        function ComputeStateForCurrentNode(this, nodeIndex)
            this.surround_nodes_index = this.surround_nodes4_each_node{nodeIndex};
            if length(this.surround_nodes_index) > this.max_ringNodes
%                 warning("length(surrounding_nodes) > max_ringNodes");
            end
            
            % 改变起始ring node，先不改变起始ring node
%             if ~this.flag_second_neighbor
%             bit = randi([1,length(surround_nodes_index)], 1, 1);
%             surrounding_node_index = circshift(surround_nodes_index, bit);
%             end
            
            this.surround_nodes_coord = this.Coord(this.surround_nodes_index, :);
            [this.State, this.RANGE_polygon] = InitializeState(this.surround_nodes_coord);  
            
           this.ring_nodes_index = this.ring_nodes4_each_node{this.nodeIndex};
           this.ring_nodes_coord = this.Coord(this.ring_nodes_index, :);
            
            if(~this.validation_flag && this.plotOptProcess)
                plot(this.Coord(nodeIndex,1), this.Coord(nodeIndex,2), 'mx')
                Plot_State_Action(this.State, [], this.RANGE_domain, this.RANGE_polygon)
            end
            
            this.State = TrimState(this.State, this.max_ringNodes);
            [~, this.neighborCells] = ComputeNodeNeighborCells(this.cellNodeTopology, this.nodeIndex);
        end
        
        function NodeIndexIncrement(this)
            this.nodeIndex = this.nodeIndex + 1;
            while sum(this.nodeIndex==this.BC_nodes) >= 1
                this.nodeIndex = this.nodeIndex + 1;
            end
        end
        
        function [observation,reward,isdone,loggedSignals] = step(this, action)
           %% 初始化值
            observation = this.State;
            reward = 0;
            isdone = false;
            loggedSignals.State = observation;
            
           %%
            if(this.validation_flag)
                if(abs(action(1))<1e-5 && abs(action(2))<1e-5)
                    disp(["action = ", num2str(action(1)), num2str(action(2))]);
                end
            end
            
            %% 当前状态对应的动作            
%             centroid = Centroid(this.ringNodesCoord);
%             centroid = Normalize(centroid, this.RANGE_polygon);
            
            valid_State = DeleteInvalidPoints(this.State);
            centroid = Centroid(valid_State);
            
            %% 如果输入的State维度大于神经网络输入维度，则采用Laplacian光顺
            if size(this.surround_nodes_coord,1) > this.max_ringNodes
                [tmp_state, ~] = InitializeState(this.surround_nodes_coord);
                centroid = Centroid(tmp_state);
                action = [0;0];
                this.lapCount = this.lapCount + 1;
%                 disp(['Using Laplacian smoothing in DRL smoothing! count = ', num2str(this.lapCount)])
            end
            
%             this.actionStorage(this.nodeIndex,:) = abs(action') ./ centroid * 100;
%             this.actionStorage(this.nodeIndex,:) = action';

            action = action + centroid';   
            
            new_point = AntiNormalize(action, this.RANGE_polygon);
            this.Coord(this.nodeIndex,:) = new_point;
            
            if(~this.validation_flag && this.plotOptProcess)
                valid_State = DeleteInvalidPoints(this.State);
                PlotCenter(valid_State, this.RANGE_polygon)
                Plot_State_Action(valid_State, action, this.RANGE_domain, this.RANGE_polygon);
            end            
            
            %% 计算当前动作的奖励
             [reward, isdone] = ComputeReward(this, new_point);
             
            if isdone 
                return;
            end
            
            NodeIndexIncrement(this);
            
           %% 计算下一步的状态
            if(this.nodeIndex>this.nNodes)
                isdone = true;
                this.IsDone = isdone;
                reward = 0;
                return;
            end
                        
            ComputeStateForCurrentNode(this, this.nodeIndex); 

            observation = this.State;
            isdone = false;
            loggedSignals.State = observation;
            this.IsDone = isdone;
        end
        
        function [reward, isdone] = ComputeReward(this, new_point)
           %% 计算动作对环境造成的改变，并评估给出奖励   
           [~, RANGE_polygon_ring] = InitializeState(this.ring_nodes_coord);  
           centroid = Centroid(this.ring_nodes_coord);
           
           vec = new_point - centroid;
           dist = sqrt(vec*vec');

            flag = SinglePointOutOfDomain(new_point, RANGE_polygon_ring, this.ring_nodes_coord);
            if flag
                isdone = true;
                this.IsDone = isdone;
                reward = this.Penalty * dist * 0.1;
                return;
            end
            
            lamda = this.lamda;            
            coeff = this.reward_coeff;          
%             [minAreaRatio, ~, averAreaRatio, minQuality, ~, averQuality] = GridQualitySummary(this.cellNodeTopology, this.Coord(:,1),this.Coord(:,2), this.Grid_stack);                        
%             reward1 = minQuality * coeff + averQuality * (1-coeff);
%             reward2 = minAreaRatio * coeff + averAreaRatio * (1-coeff);
%             reward = 0.5 *(reward1+reward2);
            
            [ShapeQuality, AngleQuality, Skewness] = GridQualityOfCells(this.neighborCells, this.Coord(:,1),this.Coord(:,2));
            minQ = ShapeQuality(1);     averQ = ShapeQuality(3);
            minAngle = AngleQuality(1); averMinAngle = AngleQuality(3);
            maxAngle = AngleQuality(2); averMaxAngle = AngleQuality(4);
            maxSkewness = Skewness(1);  averSkewness = Skewness(2);        
            % shape quality
            reward1 = coeff * minQ + (1-coeff)*averQ;
            % minimum included angle
            reward_min  = 1.0 - (60 - minAngle) / 60;                            %转换到0-1，0-最佳，1-最差
            reward_aver = 1.0 - (60 - averMinAngle) / 60;
            reward2 = ( coeff * reward_min + (1-coeff)*reward_aver );   %质量越佳，reward越大   
            % maximum included angle
            reward_max = 1.0 - (maxAngle - 60)/(180 - 60);                       %转换到0-1，0-最佳，1-最差
            reward_aver= 1.0 - (averMaxAngle - 60)/(180 - 60);          
            reward3 = ( coeff * reward_max + (1-coeff)*reward_aver );   %质量越佳，reward越大    
            % equiangle skewness， 0-最佳，1-最差
            reward_max  = 1.0 - maxSkewness;
            reward_aver = 1.0 - averSkewness;           
            reward4 = ( coeff * reward_max + (1-coeff)*reward_aver );       
            
%             reward = reward1 + reward4; % + reward2 + reward3;% 
%             reward = reward / 2;

            reward = lamda * reward1 + (1.0 - lamda)* reward4;
            
            if(~this.validation_flag && this.plotOptProcess)
                PLOT(this.Grid_stack, this.Coord(:,1), this.Coord(:,2));
            end
            
           isdone = false;
        end
        
        function ValidationCreateEnv(this, reward_coeff, validation_flag)
            this.reward_coeff = reward_coeff;
            this.validation_flag = validation_flag;
            
            PlotGrid(this, this.Coord0);
            disp('===============    原始网格    ===============');
            GridQualitySummary(this, this.Coord0);
        end
        
        function ValidationPostprocess(this)
%             figure;
%             PlotGrid(this, this.Coord);
            disp('===============    优化网格    ===============');
            GridQualitySummary(this, this.Coord);
            if(this.lapCount>0)
                disp(['DRL-Smoothing中采用Laplacian优化的次数：', num2str(this.lapCount)]);
            end
        end 
        
        function SmoothMultipleTimes(this)
            this.Coord0 = this.Coord;
        end
        
        function PlotGrid(this, Coord)
            PLOT(this.Grid_stack, Coord(:,1), Coord(:,2));            
        end
        
        function GridQualitySummary(this, Coord)
            [AreaRatioQuality, ShapeQuality, AngleQuality, Skewness] = GridQualitySummary(this.cellNodeTopology, Coord(:,1), Coord(:,2), this.Grid_stack);
            %%
            PrintGridQualitySummary(AreaRatioQuality, ShapeQuality, AngleQuality, Skewness);
        end
        
        function WriteAction2TecFile(this, actionFile)
            this.WriteTecplotTitle(actionFile);
            %% 写节点坐标
            dlmwrite(actionFile, this.Coord(:,1),'-append','delimiter','\t','precision','%.6f');
            dlmwrite(actionFile, this.Coord(:,2),'-append','delimiter','\t','precision','%.6f');
            
            %% 写节点变量数据
            dlmwrite(actionFile, this.actionStorage(:,1),'-append','delimiter','\t','precision','%.6f');
            dlmwrite(actionFile, this.actionStorage(:,2),'-append','delimiter','\t','precision','%.6f');         
            
            tmp_act = sqrt(this.actionStorage(:,1).^2 + this.actionStorage(:,2).^2);
            dlmwrite(actionFile, tmp_act,'-append','delimiter','\t','precision','%.6f');         
            
            %% 写网格拓扑
            this.WriteTecplotTopology(actionFile)
        end
        
        function WriteTecplotTitle(this, actionFile)
            title = 'Title = action data generated by Matlab';
            dlmwrite(actionFile, title, 'delimiter', '');
            
            title = 'Variables = "x", "y", "action_x", "action_y" "action"';
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            
            title = 'Zone T = "Zone1"';
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            
            title = 'ZoneType = FEPolygon';
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            
            title = ['Nodes    =', num2str(this.nNodes)];
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            
            nFaces =size(this.Grid_stack,1);
            title = ['Faces    =', num2str(nFaces)];
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            
            nCells = this.gridDataObj.nCells;
            title = ['Elements    =', num2str(nCells)];
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            
            title = ['TotalNumFaceNodes    =', num2str(2*nFaces)];
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            
            title = ['NumConnectedBoundaryFaces    =', num2str(0)];
            dlmwrite(actionFile, title, 'delimiter', '', '-append');
            title = ['TotalNumBoundaryConnections    =', num2str(0)];
            dlmwrite(actionFile, title, 'delimiter', '', '-append');            
        end
        
        function WriteTecplotTopology(this, actionFile)
            % 写face node
            dlmwrite(actionFile, this.Grid_stack(:,1:2),'-append','delimiter',' ','precision','%d');         
            % 写face left cell
            dlmwrite(actionFile, this.Grid_stack(:,3),'-append','delimiter','\t','precision','%d');         
            % 写face right cell
            dlmwrite(actionFile, this.Grid_stack(:,4),'-append','delimiter','\t','precision','%d');                              
        end
    end
end