clc;clear;
% grid_file = '1/first';
% grid_file = '5/fifth';
% grid_file = 'cylinder_m/cylinder_medium_pert';
% grid_file = 'naca0012_c/naca0012_coarse_pert';
% grid_file = '30p30n/30p30n_small_pert';
grid_file = 'rae2822/rae2822_pert';
%%
grid_path = 'D:/Codes/Mesh_Generation/meshimprove_VS2019/example/';
opt_exe = '.\tutorial-L1.exe';
chk_exe = '.\Quality.exe';
smo_exe = '.\DRL_Smoothing_1228.exe';
smo_path= 'D:/Codes/Mesh_Generation/meshimprove_VS2019/build/';
optimization_based = 0;
%%
if ~optimization_based
    pert_grid = [grid_path, grid_file, '.stl'];     %接受vtk，stl
    for method=1:5
        disp(['smooth_method = ', num2str(method)]);
        for iteration=1:10
            cmd_smo = [smo_exe, ' ',num2str(method), ' ', num2str(iteration), ' ', pert_grid, ' ', smo_path];
            status1 = system(cmd_smo, '-echo');
        end
    end
else
    %%
    pert_grid = [grid_path, grid_file, '.vtk'];     %只接受vtk
    optm_grid = [grid_path, grid_file, '_opt_smoothed.vtk'];
    for iters=1:10
        cmd_opt = [opt_exe, ' ', num2str(iters), ' ', pert_grid];
        status1 = system(cmd_opt, '-echo');
        
        if status1 == 0
            cmd_check = [chk_exe, ' ', optm_grid];
            status2 = system(cmd_check, '-echo');
        end
    end
end