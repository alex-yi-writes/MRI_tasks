%% STARTER SCRIPT : blue light task
%  Make sure you run each cell individually to avoid MATLAB throwing
%  errors.
%  In emergency contact Alex: 
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   01_04_2019      created the script

%% experiment preparation
%  always run this cell before starting an experiment.

clear all; close all; clc;
rand('state',sum(100*clock));

taskpath       = '\\fs-md\users\yiy\Documents\DOPE\BlueLight\';
savepath_pupil = [taskpath 'data\pupil\'];
savepath_behav = [taskpath 'data\behav\'];
addpath(genpath(taskpath))
addpath(genpath('\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33'))
cd(savepath_pupil)

% enter experiment details
input_prompt = {'ID (e.g. 4001)';'GROUP (1 = YAs, 2= OAs)';'1=main task, 2=immediate, 3=delayed'; ...
    'EYETRACKING? (0 = No, 1 = Yes)'; 'MRI? (0 = No, 1 = Yes)'};
defaults     = {'9999','1','1','0','0'};
input_answer = inputdlg(input_prompt, 'Input experiment details. NO SPACES!', 1, defaults);

ID           = str2num(input_answer{1,1});
AgeGroup     = str2num(input_answer{2,1});
TaskType     = str2num(input_answer{3,1});
EyeFlag      = str2num(input_answer{4,1});
MRIFlag      = str2num(input_answer{5,1});
path_ECHO    = taskpath;

%% Start Experiment
cnt = 0; errcnt = 0;
while cnt == errcnt
    try
        sca; % close all previous windows
        
        if TaskType == 1 % task
            eyetracker_filename = ['be' num2str(ID) '.edf'];
            
            [dat] = runBlueLightTask_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,path_ECHO);
            
            % save results
            save([savepath_behav num2str(ID) '_BL.mat'],'dat')
            
        elseif TaskType == 2 % tests
            eyetracker_filename = ['bi' num2str(ID) '.edf'];
            
            [dat] = runBlueLightTest_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,TaskType,path_ECHO);
            
            % save results
            save([savepath_behav num2str(ID) '_BL.mat'],'dat')
            
        elseif TaskType == 3
            eyetracker_filename = ['bd' num2str(ID) '.edf'];
            
            [dat] = runBlueLightTest_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,TaskType,path_ECHO);
            
            % save results
            save([savepath_behav num2str(ID) '_BL.mat'],'dat')
            
        else
            error('Something''s wrong...')
            
        end
    catch WhatWasTheError
        warning('error occured: trying to override')
        errcnt = errcnt+1;
    end
    cnt = cnt+1;
end

cd(savepath_pupil)
if EyeFlag ==1
    Eyelink('Message',num2str([0 200]))
    Eyelink('CloseFile'); % if interrupted, just run this and the following line manually
    Eyelink('ReceiveFile', eyetracker_filename);
    try
        fprintf('Receiving data file ''%s''\n', eyetracker_filename );
        status=Eyelink('ReceiveFile', eyetracker_filename);
        if status > 0
            fprintf('ReceiveFile status %d\n', status);
        end
        if 2==exist(dat.eyetracking.fname, 'file')
            fprintf('Data file ''%s'' can be found in ''%s''\n', eyetracker_filename, pwd );
        end
    catch
        fprintf('Problem receiving data file ''%s''\n', eyetracker_filename );
    end
    Eyelink('Shutdown')
end


%% troubleshoot

% when pupil data was not received:
% Eyelink('Initialize')
% eyetracker_filename = eyetracker_filename; % put your filename here
% Eyelink('ReceiveFile', eyetracker_filename);
