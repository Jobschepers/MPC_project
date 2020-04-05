%%%%%% Model predictive control 2020 %%%%%
%%% J.M. Schepers & M. van der Leest %%%

clc
clear all
%%

x = input('Select which matlab file to run? \n\n 1) Reference MPC \n 2) Output observer MPC \n 3) Disturbance Output MPC \n\n \n'); 

switch x
    case 1
        close all;
        run('UAV_reference.m');
        disp('Run the script again and choose option 1,2,3 or 4');
    case 2
        close all;
        run('UAV_output_dobserver.m');
        disp('Run the script again and choose option 1,2,3 or 4');
    case 3
        close all;
        run('UAV_output_distrejection.m');
        disp('Run the script again and choose option 1,2,3 or 4');
    case 4
        close all;
        run('stability_analysis.m');
        disp('Run the script again and choose option 1,2,3 or 4');
end
        