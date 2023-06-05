classdef GridDataClass < handle
    properties
        nNodes;
        nFaces;
        nCells;
        Grid_stack;
        Grid_coord;
        
        nBcTypes;
        n_of_type;
        
        BC_stack;
        BC_coord;
        BC_nodes;
        
        FAR_stack;
        FAR_coord;
        FAR_nodes;
        nFarNodes;
        
        nBoundaryNodes;
        nBoundaryFaces;
        
        Wall_stack;
        Wall_Coord;
        Wall_nodes;
        nWallNodes;
        
        cellNodeTopology;
        Ring_nodes;
        secondNeighbor;
        
        NeighborCells;
        
        STL_TR;
        delaunayTriMesh;
        invalidCellIndex;
        
        BC_WALL = 3;        %BC=3为物面
        BC_FAR  = 9;        %BC=9为外场边界
        BC_INTERIOR = 2;    %BC=2为内部面
    end
    
    methods
        function this = GridDataClass(boundaryFile, gridType)
            tmp = split(boundaryFile,'.');
            fileType = lower(tmp(end));
            if(fileType=="cas")
                [this.Grid_stack, this.Grid_coord]  = read_grid(boundaryFile, gridType);
            elseif(fileType=="stl")
                [this.Grid_stack, this.Grid_coord, this.cellNodeTopology, this.STL_TR]  = read_STL(boundaryFile);
            end
            
            this.ConstructBoundaryInfo();
        end
        
        function ConstructBoundaryInfo(this)
            this.nBcTypes = 0;
            this.n_of_type = zeros(100,1);
            this.nNodes = size(this.Grid_coord, 1);
            %%
            this.nFaces = size(this.Grid_stack, 1);
            for i = 1:this.nFaces
                bcType = this.Grid_stack(i, 7);
                if bcType ~= this.BC_INTERIOR
                    this.n_of_type(bcType) = this.n_of_type(bcType) + 1;
                end
            end
            this.nBcTypes = sum(this.n_of_type > 0);
            
            %% All BCs
            this.nBoundaryFaces = sum(this.n_of_type);
            this.BC_stack = zeros(this.nBoundaryFaces, 7);
            count = 1;
            for i = 1:this.nFaces
                bcType = this.Grid_stack(i, 7);
                if bcType ~= this.BC_INTERIOR
                    this.BC_stack(count,:) = this.Grid_stack(i, :);
                    this.BC_stack(count,3) = -1;
                    this.BC_stack(count,4) = 0;
                    
                    node1 = this.Grid_stack(i, 1);
                    node2 = this.Grid_stack(i, 2);
                    this.BC_nodes(end+1:end+2) = [node1, node2];
                    count = count + 1;
                end
            end
            this.BC_nodes = unique(this.BC_nodes);
            this.BC_coord = this.Grid_coord(this.BC_nodes, :);
            this.nBoundaryNodes = length(this.BC_nodes);
            
            %% BC_FAR
            nFarfieldFaces = this.n_of_type(this.BC_FAR);
            this.FAR_stack = zeros(nFarfieldFaces, 7);
            count = 1;
            for i = 1:this.nFaces
                bcType = this.Grid_stack(i, 7);
                if bcType == this.BC_FAR
                    this.FAR_stack(count,:) = this.Grid_stack(i, :);
                    this.FAR_stack(count,3) = -1;
                    this.FAR_stack(count,4) = 0;
                    
                    node1 = this.Grid_stack(i, 1);
                    node2 = this.Grid_stack(i, 2);
                    this.FAR_nodes(end+1:end+2) = [node1, node2];
                    count = count + 1;
                end
            end
            this.FAR_nodes = unique(this.FAR_nodes);
            this.FAR_coord = this.Grid_coord(this.FAR_nodes, :);
            this.nFarNodes = length(this.FAR_nodes);
            
            %% BC_WALL
            nWallFaces = this.n_of_type(this.BC_WALL);
            this.Wall_stack = zeros(nWallFaces, 7);
            count = 1;
            for i = 1:this.nFaces
                bcType = this.Grid_stack(i, 7);
                if bcType == this.BC_WALL
                    this.Wall_stack(count,:) = this.Grid_stack(i, :);
                    this.Wall_stack(count,3) = -1;
                    this.Wall_stack(count,4) = 0;
                    
                    node1 = this.Grid_stack(i, 1);
                    node2 = this.Grid_stack(i, 2);
                    this.Wall_nodes(end+1:end+2) = [node1, node2];
                    count = count + 1;
                end
            end
            this.Wall_nodes = unique(this.Wall_nodes);
            this.Wall_Coord = this.Grid_coord(this.Wall_nodes, :);
            this.nWallNodes = length(this.Wall_nodes);
            
            this.ConstructCell2NodeTopology()
            this.nCells = size(this.cellNodeTopology,1);
        end
        
        function initNodePositions(this, perturb_coeff)
            rng('default'); a = -1; b = 1;coeff = perturb_coeff;
            for i = 1:size(this.Grid_stack,1)
                node1 = this.Grid_stack(i,1);
                node2 = this.Grid_stack(i,2);
                len = this.Grid_stack(i,5);
                
                if sum(node1==this.BC_nodes) == 0
                    r = coeff * ( a + (b-a).*rand(1,1) );
                    this.Grid_coord(node1,1) = this.Grid_coord(node1,1) + r * len;
                    r = coeff * ( a + (b-a).*rand(1,1) );
                    this.Grid_coord(node1,2) = this.Grid_coord(node1,2) + r * len;
                end
                
                if sum(node2==this.BC_nodes) == 0
                    r = coeff * ( a + (b-a).*rand(1,1) );
                    this.Grid_coord(node2,1) = this.Grid_coord(node2,1) + r * len;
                    r = coeff * ( a + (b-a).*rand(1,1) );
                    this.Grid_coord(node2,2) = this.Grid_coord(node2,2) + r * len;
                end
            end
        end
        
        function initNodePositions_new(this, perturb_coeff)
            this.ConstructRingNodes();
            
            rng('default'); a = -1; b = 1;coeff = perturb_coeff;            
            for i=1:this.nNodes
                if sum(i==this.BC_nodes)>0
                    continue;
                end
                ringNodes = this.Ring_nodes{i};
                total = 0;
                for k = 1:length(ringNodes)
                    len = DISTANCE(i, ringNodes(k), this.Grid_coord(:,1), this.Grid_coord(:,1));
                    total = total + len;
                end
                total = total / length(ringNodes);
                
                r = coeff * ( a + (b-a).*rand(1,1) );
                this.Grid_coord(i,1) = this.Grid_coord(i,1) + r * total;
                r = coeff * ( a + (b-a).*rand(1,1) );
                this.Grid_coord(i,2) = this.Grid_coord(i,2) + r * total;                
            end
        end   
        
        function ConstructCell2NodeTopology(this)
            if isempty(this.cellNodeTopology)
%                 this.cellNodeTopology = ConstructGridTopo(this.Grid_stack);
                this.cellNodeTopology = ConstructCellTopology(this.Grid_stack);               
            end
        end
        
        function ComputeNodeNeighborCells(this)
%             nCells = size(this.cellNodeTopology,1);
%             neighborCells = cell(this.nNodes,1);
%             for nodeIndex = 1:this.nNodes
%                 neighbors = [];
%                 for i=1:nCells
%                     oneCell = this.cellNodeTopology(i,:);
%                     for j = 1:length(oneCell)
%                         node = oneCell(j);
%                         if(node == nodeIndex)
%                             neighbors(end+1,:)=oneCell;
%                             break;
%                         end
%                     end
%                 end
%                 neighborCells{nodeIndex} = cell;
%             end
        end
        
        function ConstructRingNodes(this)
            if ~isempty(this.Ring_nodes)
                return;
            end
            
            this.Ring_nodes = cell(this.nNodes, 1);
            for i=1:this.nNodes
                bcFlag = false;
                if(sum(i==this.BC_nodes)>0)
                    bcFlag = true;
                end
                
                ringNodesOfNode = NeighborNodes(i, this.Grid_stack, -1);
%                 this.Ring_nodes{i} = ringNodesOfNode;
                this.Ring_nodes{i} = SortStatePoints(this.Grid_stack, ringNodesOfNode, bcFlag);
            end
        end
        
        function ConstructNeighborNodes(this)
            if ~isempty(this.secondNeighbor)
                return;
            end
            
            this.secondNeighbor = cell(this.nNodes, 1);
            for i=1:this.nNodes               
                ringNodesOfNode = NeighborNodes(i, this.Grid_stack, -1);
                this.secondNeighbor{i} = ringNodesOfNode;
                
                if(this.ConnectedByBC(i))
                    continue;
                end
                        
                for j = 1:length(ringNodesOfNode)
                    neighborNodes = NeighborOfNeighborNodes(ringNodesOfNode(j), this.Grid_stack);
                    
                    this.secondNeighbor{i} = unique([this.secondNeighbor{i}, neighborNodes]);
                end
            end
        end
        
        function GridQualitySummaryCell(this)
            
            if(isempty(this.cellNodeTopology))
                this.ConstructCell2NodeTopology()
            end
            
            [AreaRatioQuality, ShapeQuality, AngleQuality, Skewness] = GridQualitySummary(this.cellNodeTopology, this.Grid_coord(:,1), this.Grid_coord(:,2), this.Grid_stack);
            
            PrintGridQualitySummary(AreaRatioQuality, ShapeQuality, AngleQuality, Skewness);
        end
        
        function EdgeSwap(this)
            [this.delaunayTriMesh, this.invalidCellIndex] = DelaunayMesh(this.Grid_coord(:,1), this.Grid_coord(:,2), this.Wall_nodes);
            this.Grid_coord = this.delaunayTriMesh.Points;
            count = 1;
            for i=1:size(this.delaunayTriMesh,1)
                if sum(i==this.invalidCellIndex)==0
                    this.cellNodeTopology(count,1:3) = this.delaunayTriMesh.ConnectivityList(i,:);
                    count = count + 1;
                end
            end
        end
        
        function SpringOptimize(this, iterations)
            Grid_coord_opt = SpringOptimize(this.Grid_stack, this.BC_nodes, this.Grid_coord, iterations);
            this.Grid_coord = Grid_coord_opt;
        end
        
        function SpringOptimize_1st_neighbor(this, iterations)
            this.ConstructRingNodes();
            for k=1:iterations
                disp(['Laplacian smooth iterations = ', num2str(k), ' / ', num2str(iterations)]);
                for i=1:this.nNodes
                    if(sum(this.BC_nodes==i)>0)
                        continue;
                    end

                    neigbhors = this.Ring_nodes{i};
                    neibghborCoords = this.Grid_coord(neigbhors,:);
                    
                    this.Grid_coord(i,:) = Centroid(neibghborCoords);
                end
                
                PLOT(this.Grid_stack, this.Grid_coord(:,1), this.Grid_coord(:,2), 0);
                pause(0.001);
            end
        end
        
        function SpringOptimize_2nd_neighbor(this, iterations)
            this.ConstructRingNodes();
            this.ConstructNeighborNodes();
            for k=1:iterations
                disp(['2nd neighbor Laplacian smooth iterations = ', num2str(k), ' / ', num2str(iterations)]);
                for i=1:this.nNodes
                    %对于边界点，不进行光滑
                    if(sum(this.BC_nodes==i)>0)
                        continue;
                    end

                    %对于与边界相连的点，仍然采用第一层相邻点进行光滑
                    if(this.ConnectedByBC(i))
                        neigbhors = this.Ring_nodes{i};
                    else
                        %对于内部的点，采用两层相邻点进行光滑
                        neigbhors = this.secondNeighbor{i};
                    end

                    neibghborCoords = this.Grid_coord(neigbhors,:);
                    this.Grid_coord(i,:) = Centroid(neibghborCoords);
                end

                PLOT(this.Grid_stack, this.Grid_coord(:,1), this.Grid_coord(:,2), 0);
                pause(0.001);
            end
        end
        
        function flag = ConnectedByBC(this, iNode)
            flag = false;
            for i = 1:this.nFaces
                node1 = this.Grid_stack(i,1);
                node2 = this.Grid_stack(i,2);
                
                if iNode == node1 && sum(this.BC_nodes==node2) > 0
                    flag = true;
                elseif iNode == node2 && sum(this.BC_nodes==node1) > 0
                    flag = true;
                end
            end
        end
        
        function SpringOptimizeDelaunay(this, iterations)
%             this.EdgeSwap();
%                         
%             Grid_coord_opt = SpringOptimize_delaunay(this.delaunayTriMesh, this.invalidCellIndex, this.BC_nodes, iterations, this.Grid_coord(:,1), this.Grid_coord(:,2));
%             this.Grid_coord = Grid_coord_opt;

            for i=1:iterations
                disp(['Laplacian/Edge swap smooth iterations = ', num2str(i), ' / ', num2str(iterations)]);
                this.EdgeSwap();
                Grid_coord_opt = SpringOptimize_delaunay(this.delaunayTriMesh,this.invalidCellIndex,this.BC_nodes, 1, this.Grid_coord(:,1), this.Grid_coord(:,2));
                this.Grid_coord = Grid_coord_opt;
            end
        end
        
        function GridQualitySummaryDelaunay(this)
            [AreaRatioQuality, ShapeQuality, AngleQuality, Skewness] = GridQualitySummaryDelaunay(this.delaunayTriMesh, this.invalidCellIndex, this.Grid_coord(:,1), this.Grid_coord(:,2));
            
            PrintGridQualitySummary(AreaRatioQuality, ShapeQuality, AngleQuality, Skewness);
        end
        function WriteSTL(this, filename)
            if(isempty(this.cellNodeTopology))
                this.ConstructCell2NodeTopology()
            end
            
            Points = this.Grid_coord;
            Connectivity = this.cellNodeTopology(:,1:3);
            TR = triangulation(Connectivity,Points);
            
            face_normal = TR.faceNormal;
            for i=1:size(this.cellNodeTopology,1)
                if face_normal(i,3)<0
                    tmp = this.cellNodeTopology(i,1);
                    this.cellNodeTopology(i,1) = this.cellNodeTopology(i,3);
                    this.cellNodeTopology(i,3) = tmp;
                end
            end
            Connectivity = this.cellNodeTopology(:,1:3);
            TR = triangulation(Connectivity,Points);
            
%             face_normal = TR.faceNormal;
            stlwrite(TR, filename, 'text');
            disp(['输出stl文件:', filename]);
        end
    end
end