%% STARTER SCRIPT : Scene (Indoor/Outdoor)
%  Make sure you run each cell individually to avoid MATLAB throwing
%  errors.
%  In emergency contact Alex:
%              +49 162 713 20 49
%              yeo-jin.yi.15@ucl.ac.uk
% See also runPracticeTask_MRIpilot, runMainTask_MRIpilot,
%       runMemoryTest_MRIpilot

%%  work log

%   09_10_2018      created the script

%% experiment preparation
%  always run this cell before starting an experiment.

clear all; close all; clc;
rand('state',sum(100*clock));
          
path_parent    = '\\fs-md\users\yiy\Documents\DOPE\';
taskpath       = [ path_parent 'MRPET\'];
savepath_pupil = [ path_parent 'MRPET\raw\pupil\'];
savepath_behav = [ path_parent 'MRPET\raw\behav\'];
cogentpath     = '\\fs-md\users\yiy\Documents\MATLAB\Cogent2000v1.33';
addpath(genpath(cogentpath))
% addpath(genpath('C:\Users\direkafka\Documents\MATLAB\CogGph'))
addpath(genpath(taskpath))
cd(savepath_pupil)

% enter experiment details
input_prompt = {'ID (e.g. 4001)';'GROUP (1 = YAs, 2= OAs)';'1=Practice / 2=Maintask / 3=Immediate / 4=Delayed'; ...
    'EYETRACKING? (0 = No, 1 = Yes)'; 'MRI? (0 = No, 1 = Yes)'; 'Day (1 or 2)'; 'Counterbalance Value (1 = in, 2 = out)'; 'Schedule marker (1 or 2)'};
defaults     = {'9999','1','1','0','0','0','0','0'};
input_answer = inputdlg(input_prompt, 'Input experiment details. NO SPACES!', 1, defaults);

ID           = str2num(input_answer{1,1});
AgeGroup     = str2num(input_answer{2,1});
TaskType     = str2num(input_answer{3,1});
EyeFlag      = str2num(input_answer{4,1});
MRIFlag      = str2num(input_answer{5,1});
DayNumber    = str2num(input_answer{6,1});
CnterBal     = str2num(input_answer{7,1});
Schedules    = str2num(input_answer{8,1});
PATH_echo    = path_parent;


if TaskType == 1 % practice
    input_prompt = {'Day (1 or 2)'; 'Counterbalance Value (1 = in, 2 = out)'; 'Schedule marker (1 or 2)'};
    defaults     = {'1','0','0'};
    input_answer2= inputdlg(input_prompt, 'Input experiment details for the SECOND DAY. NO SPACES!', 1, defaults);
    
    SecondDayInfo_Day            = str2num(input_answer2{1,1});
    SecondDayInfo_CounterBalance = str2num(input_answer2{2,1});
    SecondDayInfo_Schedules      = str2num(input_answer2{3,1});
end

%% Start Experiment

cnt = 0; errcnt = 0;
while cnt == errcnt
    try
        sca; % close all previous windows
        
        if TaskType == 1 % practice
                      
            eyetracker_filename = ['pr' num2str(DayNumber) num2str(ID) '.edf'];
            filestrings         = 'practice';
            
            [dat] = runPracticeTask_MRPET(ID,AgeGroup,EyeFlag,CnterBal,eyetracker_filename,DayNumber,Schedules,SecondDayInfo_CounterBalance,SecondDayInfo_Schedules,SecondDayInfo_Day,PATH_echo);
            
        elseif TaskType == 2 % main task
            eyetracker_filename = ['re' num2str(DayNumber) num2str(ID) '.edf'];
            filestrings         = 'maintask';
            load([savepath_behav num2str(ID) '_' num2str(DayNumber) '.mat']);
            
            %     if DayNumber == 1
            %         CBchar = str2num(dat.day1.RewardCategory(9));
            %     else
            %         CBchar = str2num(dat.day2.RewardCategory(9));
            %     end
            %     if CnterBal == CBchar
            [dat] = runMainTask_MRPET(ID,EyeFlag,DayNumber,MRIFlag,CnterBal,eyetracker_filename,PATH_echo,Schedules);
            %     else
            %         error('Counterbalancing information does not match! Please check again.')
            %     end
            
        elseif TaskType == 3 % immediate recall
            eyetracker_filename = ['ri' num2str(DayNumber) num2str(ID) '.edf'];
            filestrings         = 'immediate';
            load([savepath_behav num2str(ID) '_' num2str(DayNumber) '.mat']);
            
            %     if DayNumber == 1
            %         CBchar = str2num(dat.day1.RewardCategory(9));
            %     else
            %         CBchar = str2num(dat.day2.RewardCategory(9));
            %     end
            %     if CnterBal == CBchar
            %         [dat] = runMemoryTest_3Tpilot_cg_mid(ID,TaskType,EyeFlag,MRIFlag,CnterBal,eyetracker_filename,Schedules,DayNumber);
            [dat] = runMemoryTest_immediate_MRPET(ID,TaskType,EyeFlag,MRIFlag,CnterBal,eyetracker_filename,Schedules,DayNumber,PATH_echo);
            
            %     else
            %         error('Counterbalancing information does not match! Please check again.')
            %     end
            
        elseif TaskType == 4 % delayed recall
            eyetracker_filename = ['rd' num2str(DayNumber) num2str(ID) '.edf'];
            filestrings         = 'delayed';
            load([savepath_behav num2str(ID) '_' num2str(DayNumber) '.mat']);
            
            %     if DayNumber == 1
            %         CBchar = str2num(dat.day1.RewardCategory(9));
            %     else
            %         CBchar = str2num(dat.day2.RewardCategory(9));
            %     end
            %     if CnterBal == CBchar
            [dat] = runMemoryTest_delayed_MRPET(ID,TaskType,EyeFlag,MRIFlag,CnterBal,eyetracker_filename,DayNumber,Schedules,PATH_echo);
            %     else
            %         error('Counterbalancing information does not match! Please check again.')
            %     end
            
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

try
    cd(savepath_behav)
    mail = 'sfb779.behavdata@gmail.com';
    password = '*caldzn3';
    port = '465';
    setpref('Internet','SMTP_Server','smtp.gmail.com');
    setpref('Internet','E_mail',mail);
    setpref('Internet','SMTP_Username',mail);
    setpref('Internet','SMTP_Password',password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port',port);
    try
        sendmail(mail,['behavfiles_' num2str(ID) '_' num2str(DayNumber) '_' filestrings],char(datetime),[num2str(ID) '_' num2str(DayNumber) '.mat'])
    catch
        sendmail(mail,['behavfiles_' num2str(ID) '_' num2str(DayNumber) '_' filestrings],char(datestr(now)),[num2str(ID) '_' num2str(DayNumber) '.mat'])
    end
catch
    warning('data not sent! move them manually')
end


%% troubleshoot

% when pupil data was not received:
% Eyelink('Initialize')
% eyetracker_filename = eyetracker_filename; % put your filename here
% Eyelink('ReceiveFile', eyetracker_filename);
%
