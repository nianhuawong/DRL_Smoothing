classdef DRL_Opt_Action < DRL_Opt_Abstract
    
    methods
        function this = DRL_Opt_Action(max_ringNodes, boundaryFile, perturb_coeff, gridType)
            %% 动作为生成一个新点坐标
            num_out_nodes = 1;
            ActionInfo = rlNumericSpec([num_out_nodes*2 1]);%,'LowerLimit',-1, 'UpperLimit',1
            ActionInfo.Name = 'mesh DRL Action';
            ActionInfo.Description = 'new point coordinates';
            this = this@DRL_Opt_Abstract(ActionInfo, max_ringNodes, boundaryFile, perturb_coeff, gridType);
        end
    end
    
end
