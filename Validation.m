close all;
addpath(genpath('./utilities/'));
format long;
clc;
clear;
load('./agent/finalAgent_1st_20221228.mat','agent');
%%
grid_file = './grid/validation1.stl';
% grid_file = './grid/validation2.stl';
% grid_file = './grid/validation3.stl';
% grid_file = './grid/validation4.stl';
gridType      = 0;
perturb_coeff = 0.00;
reward_coeff = 0.5;
iterations   = 5;
actionFile = './action.plt';

%% 
stateDim = agent.ExperienceBuffer.ObservationDimension{1}(1);
env = DRL_Opt_Action(stateDim, grid_file, perturb_coeff, gridType);
env.ValidationCreateEnv(reward_coeff, true)
% axis([-0.25 1.25 -0.5 0.5])
% axis([-0.15 0.7 -0.25 0.25])
%%
simOptions = rlSimulationOptions('MaxSteps',env.nNodes);
if iterations > 1
    tic
    for i = 1:iterations
        disp(['iterations = ', num2str(i), ' / ', num2str(iterations)]);
        experience = sim(env,agent,simOptions);
%         env.PlotGrid(env.Coord);
        env.ValidationPostprocess();
        env.SmoothMultipleTimes();
    end
    toc
else
    experience = sim(env,agent,simOptions);
end

env.PlotGrid(env.Coord);
totalReward = sum(experience.Reward);
%%
env.ValidationPostprocess();

% env.WriteAction2TecFile(actionFile);