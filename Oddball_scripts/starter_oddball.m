%% STARTER SCRIPT : Face oddball (M/F)
%  Make sure you run each cell individually to avoid MATLAB throwing
%  errors.
%  In emergency contact Alex: 
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk

%%  work log

%   30_12_2018      created the script

%% experiment preparation
%  always run this cell before starting an experiment.

clear all; close all; clc;
rand('state',sum(100*clock));

taskpath       = '\\fs-md\users\yiy\Documents\DOPE\Oddball\';
savepath_pupil = [taskpath 'data\pupil\'];
savepath_behav = [taskpath 'data\behav\'];
addpath(genpath(taskpath))
cd(savepath_pupil)

% enter experiment details
input_prompt = {'ID (e.g. 4001)';'GROUP (1 = YAs, 2= OAs)';'1=oddball'; ...
    'EYETRACKING? (0 = No, 1 = Yes)'; 'MRI? (0 = No, 1 = Yes)'; 'Schedule marker (1 or 2)'};
defaults     = {'9999','1','1','1','1','0'};
input_answer = inputdlg(input_prompt, 'Input experiment details. NO SPACES!', 1, defaults);

ID           = str2num(input_answer{1,1});
AgeGroup     = str2num(input_answer{2,1});
TaskType     = str2num(input_answer{3,1});
EyeFlag      = str2num(input_answer{4,1});
MRIFlag      = str2num(input_answer{5,1});
Schedules    = str2num(input_answer{6,1});

%% Start Experiment

sca; % close all previous windows

if TaskType == 1 % oddball
    eyetracker_filename = ['ob' num2str(ID) '.edf'];
    
    [dat] = runOddball_cg(ID,AgeGroup,EyeFlag,MRIFlag,eyetracker_filename,Schedules);
    
    % save results
    save([savepath_behav num2str(ID) '_OB.mat'],'dat')
    
else
    error('Something''s wrong...')
    
end


cd(savepath_pupil)
if EyeFlag ==1
    Eyelink('Message',num2str([0 200]))
    Eyelink('CloseFile'); % if interrupted, just run this and the following line manually6
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
