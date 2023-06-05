clear;clc;close all;
addpath(genpath('./utilities'));
%% 环境、状态及动作定义
env = DRL_Opt_Action(8, '', -1, -1);   %当取1层邻点时，数量最大取8；当取2层相邻点时，数量最大取20
obsInfo = getObservationInfo(env); 
actInfo = getActionInfo(env);
rng('default')
%% 建立critic网络，DQN和DDPG将观察值state和动作值action同时作为Critic输入
L1 = 16;    %当取1层邻点时，取16；当取2层相邻点时，取40
statePath = [imageInputLayer([obsInfo.Dimension(1) obsInfo.Dimension(2) 1], 'Normalization', 'none', 'Name', 'State')
             fullyConnectedLayer(4*L1,'Name','CriticStateFC1')
             reluLayer('Name','CriticStateRelu1')
             fullyConnectedLayer(L1,'Name','CriticStateFC2')
             reluLayer('Name','CriticStateRelu2')
             fullyConnectedLayer(L1/2,'Name','CriticStateFC3')
%              reluLayer('Name','CriticStateRelu3')
%              fullyConnectedLayer(L1/4,'Name','CriticStateFC4')
             ];
actionPath = [imageInputLayer([actInfo.Dimension(1) actInfo.Dimension(2) 1], 'Normalization', 'none', 'Name', 'Action')
%               fullyConnectedLayer(L1,'Name','CriticActionFC1')
              reluLayer('Name','CriticActionRelu1')
              fullyConnectedLayer(L1/2,'Name','CriticActionFC2')
             ];
commonPath = [concatenationLayer(3, 2, 'Name', 'concat')
              %additionLayer(2,'Name','add')
%               reluLayer('Name','CriticCommonRelu1')
              fullyConnectedLayer(L1/4,'Name','CriticCommonFC1')
%               reluLayer('Name','CriticCommonRelu2')
              fullyConnectedLayer(1,'Name','output')
              ];
          
criticNetwork = layerGraph(statePath);
criticNetwork = addLayers(criticNetwork, actionPath);
criticNetwork = addLayers(criticNetwork, commonPath);
criticNetwork = connectLayers(criticNetwork, 'CriticStateFC3', 'concat/in1');   %'add/in1'
criticNetwork = connectLayers(criticNetwork, 'CriticActionFC2', 'concat/in2');  %'concat/in2'
% plot(criticNetwork)
criticOpts = rlRepresentationOptions('LearnRate',5e-4,'GradientThreshold', 1);%, 'Optimizer',"rmsprop"
critic = rlQValueRepresentation(criticNetwork,obsInfo,actInfo,'Observation',{'State'},'Action',{'Action'},criticOpts);

%% 建立actor网络，将观察state作为输入
actorNetwork = [
    imageInputLayer([obsInfo.Dimension(1) obsInfo.Dimension(2) 1],'Normalization','none','Name','state')
%     fullyConnectedLayer(2*L1,'Name','ActorFC1')
%     leakyReluLayer('Name','ActorRelu1')    
%     fullyConnectedLayer(L1,'Name','ActorFC2')
%     leakyReluLayer('Name','ActorRelu2')    
%     fullyConnectedLayer(L1,'Name','ActorFC3')
%     leakyReluLayer('Name','ActorRelu3')

    fullyConnectedLayer(2*L1,'Name','ActorFC1')
    reluLayer('Name','ActorRelu1')
    fullyConnectedLayer(L1,'Name','ActorFC2')
    reluLayer('Name','ActorRelu2')
    fullyConnectedLayer(L1,'Name','ActorFC3') 
    reluLayer('Name','ActorRelu3')
%     fullyConnectedLayer(16,'Name','ActorFC4')
%     reluLayer('Name','ActorRelu4')
%     fullyConnectedLayer(8,'Name','ActorFC5')
%     reluLayer('Name','ActorRelu5')

    fullyConnectedLayer(actInfo.Dimension(1),'Name','actor')
%     tanhLayer('Name','actor')
    ];
% plot(layerGraph(actorNetwork))
actorOpts = rlRepresentationOptions('LearnRate',1e-5,'GradientThreshold', 1);
actor = rlDeterministicActorRepresentation(actorNetwork,obsInfo,actInfo,'Observation',{'state'}, 'Action', {'actor'}, actorOpts);

%% 建立智能体DDPG agent
agentOpts = rlDDPGAgentOptions(...   
    'TargetSmoothFactor',1e-3,...
    'ExperienceBufferLength',1e6,...
    'DiscountFactor',0.99,...
    'MiniBatchSize',64,...
    'SampleTime', 1, ...
    'SaveExperienceBufferWithAgent', true,...
    'NumStepsToLookAhead', 1);
agentOpts.NoiseOptions.Variance = 0.3;
agentOpts.NoiseOptions.VarianceDecayRate = 1e-4;

agent = rlDDPGAgent(actor,critic,agentOpts);

%% 训练智能体
averQuality = 0.95;
steps = env.nNodes;
trainOpts = rlTrainingOptions(...
    'MaxEpisodes',50000,...
    'SaveAgentCriteria',"EpisodeCount",...
    'SaveAgentValue',10000,...
    'SaveAgentDirectory', pwd + "\agent\AutoSave",...
    'MaxStepsPerEpisode',steps,...
    'Verbose',true,...
    'Plots','none',...      %training-progress/none
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',averQuality*steps,...
    'ScoreAveragingWindowLength',10); 

%% 是否加载预训练的agent
loadAgent = 0;
if loadAgent
    load('./agent/finalAgent_1st_20221228.mat','agent');
    env.max_ringNodes = agent.ExperienceBuffer.ObservationDimension{1}(1);
end

%% 是否训练agent
doTraining = true;
if doTraining  
    trainingStats = train(agent,env,trainOpts);
    save("./agent/finalAgent_"+num2str(steps)+".mat",'agent')
end

% save("./agent/finalAgent_1st_930.mat",'agent')